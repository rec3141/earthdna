#!/usr/bin/env python3
"""FASTA -> DNABERT-2 / DNABERT-S embeddings (mean of the last hidden layer over tokens) as a fastsne .fsb file.

DNABERT-S (zhihan1996/DNABERT-S) is DNABERT-2 fine-tuned so fragments of one species embed together, the closest
off-the-shelf model to what the globes ask for; DNABERT-2 (zhihan1996/DNABERT-2-117M) is the plain language model.
Both use ALiBi, so 4 kb windows (~900 BPE tokens) fit without truncation.

Runs on one GPU with uv (no install):
  uv run --python 3.10 --with torch==2.3.1 --with transformers==4.29.2 --with einops --with numpy \
      pipeline/embed_dnabert.py windows.fa out.fsb --model zhihan1996/DNABERT-S [--batch 32] [--max-tokens 2048]
The models' bundled Triton flash attention needs an old Triton; it is switched off (plain attention is used).
"""
import argparse, gzip, os, struct, sys, time
import numpy as np, torch
from transformers import AutoModel, AutoTokenizer


def read_fasta(path):
    """Read plain or .gz FASTA into aligned (ids, seqs); reject empty records or malformed headers."""
    ids, seqs, sid, cur = [], [], None, None
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        for l in f:
            if l.startswith(">"):
                if cur is not None:
                    seq = "".join(cur)
                    if not seq:
                        raise SystemExit(f"empty sequence for {sid} in {path}")
                    ids.append(sid)
                    seqs.append(seq)
                toks = l[1:].split()
                if not toks:
                    raise SystemExit(f"empty FASTA header in {path}")
                sid, cur = toks[0], []
            else:
                s = l.strip()
                if not s:
                    continue
                if cur is None:
                    raise SystemExit(f"sequence data before first FASTA header in {path}")
                cur.append(s.upper().replace("U", "T"))
        if cur is not None:
            seq = "".join(cur)
            if not seq:
                raise SystemExit(f"empty sequence for {sid} in {path}")
            ids.append(sid)
            seqs.append(seq)
    return ids, seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta"); ap.add_argument("out")
    ap.add_argument("--model", default="zhihan1996/DNABERT-S")
    ap.add_argument("--batch", type=int, default=32); ap.add_argument("--max-tokens", type=int, default=2048)
    a = ap.parse_args()

    ids, seqs = read_fasta(a.fasta)

    tok = AutoTokenizer.from_pretrained(a.model, trust_remote_code=True)
    model = AutoModel.from_pretrained(a.model, trust_remote_code=True)
    for name, mod in sys.modules.items():          # disable the bundled Triton kernel
        if name.endswith("bert_layers") and hasattr(mod, "flash_attn_qkvpacked_func"): mod.flash_attn_qkvpacked_func = None
    model = model.cuda().eval()

    order = np.argsort([len(s) for s in seqs])     # similar lengths per batch: less padding
    out = np.zeros((len(seqs), model.config.hidden_size), np.float32)
    t0 = time.time()
    with torch.no_grad(), torch.autocast("cuda", dtype=torch.float16):
        for k in range(0, len(seqs), a.batch):
            ix = order[k:k + a.batch]
            enc = tok([seqs[i] for i in ix], return_tensors="pt", padding=True, truncation=True, max_length=a.max_tokens)
            h = model(enc["input_ids"].cuda(), attention_mask=enc["attention_mask"].cuda())[0].float()
            m = enc["attention_mask"].cuda().bool().unsqueeze(-1)
            denom = m.sum(1).clamp_min(1).float()
            out[ix] = (torch.where(m, h, torch.zeros_like(h)).sum(1) / denom).cpu().numpy()
            if k // a.batch % 50 == 0:
                print(f"{k + len(ix):,}/{len(seqs):,}  {(k + len(ix)) / max(time.time() - t0, 1e-9):.0f} seq/s", flush=True)

    assert out.shape == (len(ids), model.config.hidden_size)
    tmp_out = a.out + ".tmp"
    with open(tmp_out, "wb") as o:
        o.write(b"FSB1" + struct.pack("<iq", out.shape[1], len(ids)))
        o.write(out.tobytes())
        o.write("\n".join(ids).encode())
    os.replace(tmp_out, a.out)
    print(f"{len(ids):,} sequences -> {a.out} ({a.model})")


if __name__ == "__main__":
    main()
