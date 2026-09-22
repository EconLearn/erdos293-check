#!/usr/bin/env python3
"""
check_v9_certs.py FILE [LO HI] -- standalone check of the v(9) certificate file.
Jude Wallis, 2026-09-21.

For each line 'm: n_1 ... n_9': the n_i are 9 distinct integers >= 2, m is one of them, and
sum 1/n_i == 1, computed with fractions.Fraction. Also checks that every m in [LO, HI]
(default [2, 3000000]) has exactly one line. If all checks pass, every integer in [LO, HI]
occurs in some 9-term representation of 1, so v(9) > HI.
"""
import sys
from fractions import Fraction

fn = sys.argv[1]
lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (2, 3000000)
seen = bytearray(hi + 1)
bad = 0
for line in open(fn):
    if line.startswith("#") or not line.strip():
        continue
    head, rest = line.split(":")
    m, ns = int(head), [int(x) for x in rest.split()]
    ok = len(ns) == 9 and len(set(ns)) == 9 and min(ns) >= 2 and m in ns and sum(Fraction(1, n) for n in ns) == 1
    if not ok:
        bad += 1
        print("BAD:", line.strip())
    elif lo <= m <= hi:
        if seen[m]:
            bad += 1
            print("DUPLICATE:", m)
        seen[m] = 1
gaps = [m for m in range(lo, hi + 1) if not seen[m]]
print(f"bad lines: {bad}; integers in [{lo}, {hi}] without a certificate: {len(gaps)} {gaps[:10]}")
if bad == 0 and not gaps:
    print(f"OK: every m in [{lo}, {hi}] occurs in a 9-term representation of 1, so v(9) > {hi}")
