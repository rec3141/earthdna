#!/usr/bin/env python3
"""FASTA -> Evo 2 embeddings as fastsne .fsb tables (float32 rows in FASTA order, ids after the data).

One forward pass gives every requested layer x pooling combination, written as OUT.<layer>.<pool>.fsb:
  mean  average over all positions
  tail  average over the second half only (Goodfire's phylogeny-manifold recipe: positions with context)
--rc also writes <pool>rc variants: the average with the reverse complement's vector (Evo 2 is causal, so strand
matters; doubles cost). The plain variants are still written.

Sequences are length-sorted within chunks and right-padded into batches: the model is causal, so padding after
a sequence cannot change its activations, and pooling stops at the real length. --check compares a few batched
vectors with unbatched ones. Each finished chunk is saved under OUT.parts/ so a killed job resumes;
--start/--end select a shard of chunks for parallel jobs, and --assemble joins all shards into the .fsb files.

usage (inside the evo2 container):
  embed_evo2.py in.fasta OUT --local_path evo2_7b.pt [--layers blocks.24.mlp.l3,blocks.28.mlp.l3]
                [--pools mean,tail] [--rc] [--start C --end C] [--assemble]
"""
import argparse, gzip, os, struct, time
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("fasta"); ap.add_argument("out")
ap.add_argument("--model_name", default="evo2_7b")
ap.add_argument("--local_path")
ap.add_argument("--layers", default="blocks.24.mlp.l3,blocks.28.mlp.l3")
ap.add_argument("--pools", default="mean,tail")
ap.add_argument("--rc", action="store_true")
ap.add_argument("--chunk", type=int, default=2048, help="sequences per saved part")
ap.add_argument("--tokens", type=int, default=16384, help="padded tokens per batch")
ap.add_argument("--start", type=int, default=0, help="first chunk of this shard")
ap.add_argument("--end", type=int, default=None, help="one past the last chunk of this shard")
ap.add_argument("--check", type=int, default=0, help="compare N batched vectors with unbatched ones")
ap.add_argument("--assemble", action="store_true", help="join finished parts into the .fsb files and exit")
ap.add_argument("--no-kernels", action="store_true", help="plain vortex instead of its opt-in Triton HC kernels (as\n               evo2's own test runs; the kernels gave NaN on a model split over two GPUs)")
ap.add_argument("--no-fp8", action="store_true", help="bf16 input projections instead of FP8 (7B only; for GPUs below\n               compute capability 8.9 such as the A30; embeddings differ slightly from FP8 runs)")
a = ap.parse_args()
layers, pools = a.layers.split(","), a.pools.split(",")
combos = [(l, p) for l in layers for p in pools]
if a.rc: combos += [(l, p + "rc") for l in layers for p in pools]
parts = a.out + ".parts"
os.makedirs(parts, exist_ok=True)

ids, seqs, cur = [], [], None
with (gzip.open if a.fasta.endswith(".gz") else open)(a.fasta, "rt") as f:
    for line in f:
        if line[0] == ">":
            if cur is not None: seqs.append("".join(cur))
            ids.append(line[1:].split()[0]); cur = []
        else:
            cur.append(line.strip().upper().replace("U", "T"))
    if cur is not None: seqs.append("".join(cur))
nchunk = (len(seqs) + a.chunk - 1) // a.chunk
part = lambda c: f"{parts}/{c:06d}.npz"
tag = lambda l, p: l.replace(".", "_") + "." + p

if a.assemble:
    missing = [c for c in range(nchunk) if not os.path.exists(part(c))]
    if missing: raise SystemExit(f"{len(missing)} of {nchunk} chunks missing, first {missing[0]}")
    outs = {k: open(f"{a.out}.{tag(*k)}.fsb", "wb") for k in combos}
    for c in range(nchunk):
        z = np.load(part(c))
        for k, o in outs.items():
            if c == 0: o.write(b"FSB1" + struct.pack("<iq", z[tag(*k)].shape[1], len(ids)))
            o.write(z[tag(*k)].astype(np.float32).tobytes())
    for o in outs.values(): o.write("\n".join(ids).encode()); o.close()
    print(f"{len(ids)} sequences -> {', '.join(o.name for o in outs.values())}")
    raise SystemExit

import torch
from evo2 import Evo2
if a.no_fp8:
    import evo2.models
    evo2.models.HAS_TE = False   # evo2's own fallback: a 7B without Transformer Engine runs bf16 projections
import torch
# vortex's opt-in Triton kernels return NaN for every layer past the first GPU when a model is split over several
# (40B on 2 x H100: loss NaN with them, 0.216 = evo2's reference without), so they only run on a single GPU
model = Evo2(a.model_name, local_path=a.local_path, use_kernels=not a.no_kernels and torch.cuda.device_count() == 1)
dev = "cuda:0"
comp = str.maketrans("ACGTRYKMBVDHN", "TGCAYRMKVBHDN")

@torch.inference_mode()
def embed(batch):
    """list of str -> {(layer, pool): float32 [n, D]}, batches padded to the longest sequence"""
    L = max(len(s) for s in batch)
    x = torch.full((len(batch), L), ord("N"), dtype=torch.int)
    for i, s in enumerate(batch):
        x[i, :len(s)] = torch.tensor(model.tokenizer.tokenize(s), dtype=torch.int)
    # our own hooks, copying each output as it is produced: evo2's return_embeddings keeps a bare reference, and on a
    # model split over two GPUs the layers past the first GPU came back as all zeros (40B, 2026-11-09)
    emb, handles = {}, []
    def keep(name):
        def hook(_, __, out): emb[name] = (out[0] if isinstance(out, tuple) else out).detach().float().clone()
        return hook
    for l in layers: handles.append(model.model.get_submodule(l).register_forward_hook(keep(l)))
    try: model.model.forward(x.to(dev))
    finally:
        for hd in handles: hd.remove()
    n = torch.tensor([len(s) for s in batch], device=dev)
    pos = torch.arange(L, device=dev)[None, :]
    w = {"mean": (pos < n[:, None]), "tail": (pos < n[:, None]) & (pos >= (n // 2)[:, None])}
    # "...n" pools (meann, tailn) scale every position to unit length first: in the 40B the residual stream jumps to
    # ~1e16 at block 23 and only its direction carries information from there on
    w.update({k + "n": v for k, v in list(w.items())})
    out = {}
    for l in layers:
        h = emb[l]                  # a model split over GPUs returns each layer on its own device
        hn = None
        for p in pools:
            m = w[p].to(h.device).float()[..., None]
            if p.endswith("n"):
                if hn is None: hn = h / h.norm(dim=-1, keepdim=True).clamp_min(1e-30)
                x = hn
            else: x = h
            out[l, p] = ((x * m).sum(1) / m.sum(1)).cpu().numpy()
    return out

def safe_embed(batch):
    """embed(), halving the batch when the GPU runs out of room (cuFFT plans fail with CUFFT_INTERNAL_ERROR
    rather than an out-of-memory error when a padded batch of long sequences is too big)"""
    try:
        return embed(batch)
    except (RuntimeError, torch.OutOfMemoryError) as e:
        if len(batch) == 1 or not ("cuFFT" in str(e) or "out of memory" in str(e)): raise
        # cuFFT keeps plans (and their workspace) per device; drop them so the halves can allocate
        for dv in range(torch.cuda.device_count()):
            torch.backends.cuda.cufft_plan_cache[dv].clear()
        torch.cuda.empty_cache()
        print(f"batch of {len(batch)} x {max(map(len, batch))} nt failed ({str(e)[:40]}); halving", flush=True)
        h = len(batch) // 2
        x, y = safe_embed(batch[:h]), safe_embed(batch[h:])
        return {k: np.concatenate([x[k], y[k]]) for k in x}

def run(chunk):
    order = sorted(range(len(chunk)), key=lambda i: len(chunk[i]))
    res = {k: np.zeros((len(chunk), 0), np.float32) for k in combos}
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and (j + 1 - i) * len(chunk[order[j]]) <= a.tokens: j += 1
        idx = order[i:j]
        e = safe_embed([chunk[t] for t in idx])
        if a.rc:
            r = safe_embed([chunk[t].translate(comp)[::-1] for t in idx])
            e.update({(l, p + "rc"): (e[l, p] + r[l, p]) / 2 for l, p in list(e)})
        for k, v in e.items():
            if res[k].shape[1] == 0: res[k] = np.zeros((len(chunk), v.shape[1]), np.float32)
            res[k][idx] = v
        i = j
    return res

if a.check:
    test = seqs[:a.check]
    b = run(test)
    for i, s in enumerate(test):
        u = embed([s])
        for k in u:
            cos = float(np.dot(b[k][i], u[k][0]) / np.linalg.norm(b[k][i]) / np.linalg.norm(u[k][0]))
            if cos < 0.999: print(f"WARN batched vs unbatched cosine {cos:.5f} for {ids[i]} {k}")
    print(f"check done on {len(test)} sequences", flush=True)

end = min(a.end if a.end is not None else nchunk, nchunk)
todo = [c for c in range(a.start, end) if not os.path.exists(part(c))]
print(f"{len(seqs)} sequences, chunks {a.start}..{end - 1}, {len(todo)} to do", flush=True)
t0, done = time.time(), 0
for c in todo:
    chunk = seqs[c * a.chunk:(c + 1) * a.chunk]
    res = run(chunk)
    np.savez(part(c) + ".tmp.npz", **{tag(*k): v.astype(np.float32) for k, v in res.items()})   # float32: 40B residual channels reach ~1e16
    os.replace(part(c) + ".tmp.npz", part(c))
    done += len(chunk)
    rate = done / (time.time() - t0)
    left = sum(len(seqs[x * a.chunk:(x + 1) * a.chunk]) for x in todo) - done
    print(f"chunk {c}  {rate:.1f} seq/s  eta {left / rate / 3600:.1f} h", flush=True)
