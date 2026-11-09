#!/usr/bin/env python3
"""Score how much a genome contig raises Evo 2's likelihood of each candidate 16S (ssu_link_set.py output).

For every (test genome, candidate 16S) pair:
  cond    sum of log p(16S token | contig + preceding 16S tokens), over 16S tokens 2..L
  uncond  sum of log p(16S token | preceding 16S tokens) for the 16S alone, same tokens
  gain    cond - uncond: what the contig adds (a pointwise mutual information); the genome's own 16S should gain most
Both orientations of the contig are scored (contig strand relative to the rRNA operon is unknown); gain takes the
larger. Writes OUT.tsv: test, candidate, relation, len, uncond, cond_fwd, cond_rev, gain.

usage (in the evo2 container): ssu_link_score.py SET OUT --model_name evo2_40b --local_path evo2_40b.pt
  SET = ssu_link_set.py prefix (SET.contigs.fa, SET.ssu.fa, SET.pairs.tsv)
"""
import argparse, time
import torch
from evo2 import Evo2

ap = argparse.ArgumentParser()
ap.add_argument("set"); ap.add_argument("out")
ap.add_argument("--model_name", default="evo2_40b"); ap.add_argument("--local_path", required=True)
ap.add_argument("--start", type=int, default=0); ap.add_argument("--end", type=int, default=10 ** 9)
a = ap.parse_args()

def fasta(p):
    d, n = {}, None
    for l in open(p):
        if l[0] == ">": n = l[1:].strip(); d[n] = []
        else: d[n].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}
contig, ssu = fasta(a.set + ".contigs.fa"), fasta(a.set + ".ssu.fa")
pairs = [l.rstrip("\n").split("\t") for l in open(a.set + ".pairs.tsv")][1:]
tests = sorted({p[0] for p in pairs})[a.start:a.end]
rc = lambda s: s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

model = Evo2(a.model_name, local_path=a.local_path)
dev = next(model.model.parameters()).device

@torch.no_grad()
def tail_logprob(prefix, target):
    """sum log p(target[1:] | prefix + target[:i]) for i = 1..len-1"""
    ids = torch.tensor(model.tokenizer.tokenize(prefix + target), dtype=torch.int).unsqueeze(0).to(dev)
    logits, _ = model.model.forward(ids)            # as evo2's own test does; logits sit on the last GPU
    lp = torch.log_softmax(logits[0, :-1].float(), -1)
    nxt = ids[0, 1:].long().to(lp.device)
    tok = lp.gather(1, nxt.unsqueeze(1)).squeeze(1)
    return float(tok[len(prefix):].sum())   # tok[i] predicts token i+1: target[1:] is predicted from index len(prefix)

uncond = {}
out = open(a.out + ".tsv", "w")
out.write("test\tcandidate\trelation\tlen\tuncond\tcond_fwd\tcond_rev\tgain\n")
t0 = time.time()
for k, g in enumerate(tests):
    c = contig[g]
    for p in (p for p in pairs if p[0] == g):
        s = ssu[p[1]]
        if p[1] not in uncond: uncond[p[1]] = tail_logprob("", s)
        f, r = tail_logprob(c, s), tail_logprob(rc(c), s)
        out.write(f"{g}\t{p[1]}\t{p[2]}\t{len(s)}\t{uncond[p[1]]:.3f}\t{f:.3f}\t{r:.3f}\t{max(f, r) - uncond[p[1]]:.3f}\n")
    out.flush()
    print(f"{k + 1}/{len(tests)} genomes, {time.time() - t0:.0f} s", flush=True)
