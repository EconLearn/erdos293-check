#!/usr/bin/env python3
"""
split_stage.py -- certificates by splitting one term.  Jude Wallis, 2026-09-21.

If d < m, s = m - d divides d^2 (equivalently s | m^2), and d occurs in an 8-term representation W
of 1, then 1/d = 1/m + 1/e with e = d*m/s, so W - {d} + {m, e} is a 9-term representation
containing m, provided m and e are not already in W (m never is when m is missing at k = 8).

  split_stage.py plan  FAILS BITMAP OPTIONS  > plan.txt      (m d_1 d_2 ... : up to OPTIONS choices of d,
                                                              d present at k = 8, largest d first)
  split_stage.py build plan.txt W8FILE... > certs.txt        (W8FILE: findm 8 output for the d's)
"""
import sys
from math import lcm


def spf_table(L):
    spf = list(range(L + 1))
    for i in range(2, int(L ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, L + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def divisors_of_square(m, spf):
    f = {}
    while m > 1:
        p = spf[m]; f[p] = f.get(p, 0) + 1; m //= p
    ds = [1]
    for p, e in f.items():
        ds = [d * p ** k for d in ds for k in range(2 * e + 1)]
    return ds


if sys.argv[1] == "plan":
    fails = [int(x) for x in open(sys.argv[2]).read().split()]
    bm = open(sys.argv[3], "rb").read()
    seen = lambda n: (bm[n >> 3] >> (n & 7)) & 1
    k = int(sys.argv[4])
    spf = spf_table(max(fails) + 1)
    for m in fails:
        ds = sorted((m - s for s in divisors_of_square(m, spf) if s < m - 1 and 2 * s != m and seen(m - s)), reverse=True)
        print(m, *ds[:k])

elif sys.argv[1] == "build":
    plan = [list(map(int, l.split())) for l in open(sys.argv[2]) if l.strip()]
    W = {}
    for fn in sys.argv[3:]:
        for line in open(fn):
            head, _, rest = line.partition(":")
            if "FAIL" not in rest and rest.strip():
                W[int(head)] = [int(x) for x in rest.split()]
    made = missed = 0
    for row in plan:
        m, opts = row[0], row[1:]
        for d in opts:
            w = W.get(d)
            if w is None:
                continue
            s = m - d
            e = d * m // s
            if d * m % s or m in w or e in w:
                continue
            c = sorted([x for x in w if x != d] + [m, e])
            L = lcm(*c)
            assert len(set(c)) == 9 and sum(L // x for x in c) == L
            print(f"{m}: " + " ".join(map(str, c)))
            made += 1
            break
        else:
            missed += 1
    print(f"made {made}, missed {missed}", file=sys.stderr)
