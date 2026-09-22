#!/usr/bin/env python3
"""padic_build.py P JOBFILE OUTFILE [JOBFILE OUTFILE ...] -- turn findfrac hits into 9-term certificates for P.
Jude Wallis, 2026-09-22. Each certificate is re-added with fractions.Fraction and marked OK or BAD."""
import sys
from fractions import Fraction
P = int(sys.argv[1]); pairs = sys.argv[2:]
certs = []
for jf, of in zip(pairs[::2], pairs[1::2]):
    for job, out in zip(open(jf), open(of)):
        if "FAIL" in out: continue
        f = list(map(int, job.split()))
        extra = f[5:5 + f[3] - 1]                 # forbidden list minus P = the other multiples of P
        u = list(map(int, out.split(":")[1].split()))
        c = sorted([P] + extra + u)
        ok = len(c) == 9 and len(set(c)) == 9 and sum(Fraction(1, n) for n in c) == 1
        certs.append((ok, c))
for ok, c in certs:
    print(("OK  " if ok else "BAD ") + f"{P}: " + " ".join(map(str, c)))
