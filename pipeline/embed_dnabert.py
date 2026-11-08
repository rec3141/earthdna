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
import argparse, struct, sys, time
import numpy as np, torch
from transformers import AutoModel, AutoTokenizer

ap = argparse.ArgumentParser()
ap.add_argument("fasta"); ap.add_argument("out")
ap.add_argument("--model", default="zhihan1996/DNABERT-S")
ap.add_argument("--batch", type=int, default=32); ap.add_argument("--max-tokens", type=int, default=2048)
a = ap.parse_args()

ids, seqs, cur = [], [], []
for l in open(a.fasta):
    if l[0] == ">":
        if cur: seqs.append("".join(cur))
        ids.append(l[1:].split()[0]); cur = []
    else: cur.append(l.strip().upper())
if cur: seqs.append("".join(cur))

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
        m = enc["attention_mask"].cuda().unsqueeze(-1).float()
        out[ix] = ((h * m).sum(1) / m.sum(1)).cpu().numpy()
        if k // a.batch % 50 == 0:
            print(f"{k + len(ix):,}/{len(seqs):,}  {(k + len(ix)) / (time.time() - t0):.0f} seq/s", flush=True)
with open(a.out, "wb") as o:
    o.write(b"FSB1" + struct.pack("<iq", out.shape[1], len(ids))); o.write(out.tobytes()); o.write("\n".join(ids).encode())
print(f"{len(ids):,} sequences -> {a.out} ({a.model})")
