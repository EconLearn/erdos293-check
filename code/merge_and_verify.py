#!/usr/bin/env python3
"""
merge_and_verify.py -- merge the per-worker outputs of ef_enum and check them.
Jude Wallis, 2026-09-19. Written with the help of an AI coding assistant; checked and run by me.

  merge_and_verify.py PREFIX K NWORKERS [EXPECTED_COUNT]
  merge_and_verify.py --parts K EXPECTED_COUNT OUTPREFIX prefixA:w,w,... prefixB:w,...

Checks performed (all exact):
  1. total number of representations == EXPECTED_COUNT (A006585(K)) if given;
  2. v(K) = least m >= 2 whose bit is not set in the OR of all worker bitmaps;
  3. every m in 2 .. v(K)-1 has an explicit witness; each witness is re-verified from
     scratch with Python integers: K distinct positive integers, contains m, and the
     unit fractions sum to exactly 1 (common-denominator integer identity, no floats);
  4. no witness anywhere contains v(K) (sanity).
Writes PREFIX.kK.witnesses.txt (one verified witness per m < v(K)).
"""
import sys
from fractions import Fraction

if sys.argv[1] == "--parts":
    # merge_and_verify.py --parts K EXPECTED OUTPREFIX prefixA:w,w,w prefixB:w,w ...
    # (used when one slow worker's share was re-run as several finer shards; the shards
    #  idx % 81 in {0, 9, ..., 72} are exactly the prefixes idx % 9 == 0)
    K, expected, prefix = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    bases = [f"{spec.split(':')[0]}.k{K}.w{w}" for spec in sys.argv[5:] for w in spec.split(':')[1].split(',')]
else:
    prefix, K, NW = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    expected = int(sys.argv[4]) if len(sys.argv) > 4 else None
    bases = [f"{prefix}.k{K}.w{w}" for w in range(NW)]

total, maxq = 0, 0
bitmap = None
wit = {}
for base in bases:
    lines = open(base + ".count").read().split()
    total += int(lines[0]); maxq = max(maxq, int(lines[2]))
    b = open(base + ".bitmap", "rb").read()
    bitmap = bytearray(b) if bitmap is None else bytearray(x | y for x, y in zip(bitmap, b))
    for line in open(base + ".witness"):
        m, rest = line.split(":")
        m = int(m)
        if m not in wit:
            wit[m] = [int(x) for x in rest.split()]

print(f"k={K}: total representations = {total}" + (f"  (expected {expected}: {'MATCH' if total == expected else 'MISMATCH'})" if expected is not None else ""))
print(f"largest remaining-fraction denominator q seen: {maxq}")

def seen(m):
    return (bitmap[m >> 3] >> (m & 7)) & 1

v = 2
while seen(v):
    v += 1
print(f"v({K}) = {v}   (least integer >= 2 whose bit is unset; bitmap covers m < {len(bitmap) * 8})")

# gaps just above v(K), for the record
nxt = [m for m in range(v, v + 200000) if m < len(bitmap) * 8 and not seen(m)][:10]
print(f"first absent integers >= 2: {nxt}")

bad = 0
for m in range(2, v):
    s = wit.get(m)
    ok = (s is not None and len(s) == K and len(set(s)) == K and m in s and min(s) >= 2
          and sum(Fraction(1, n) for n in s) == 1)
    if not ok:
        bad += 1
        print("BAD/MISSING witness for", m, s)
print(f"witnesses verified exactly for m = 2..{v - 1}: {v - 2 - bad} ok, {bad} bad")
contains_v = [m for m, s in wit.items() if v in s]
print(f"witnesses containing {v}: {len(contains_v)} (must be 0)")

with open(f"{prefix}.k{K}.witnesses.txt", "w") as f:
    f.write(f"# One explicit representation of 1 as a sum of {K} distinct unit fractions containing m, for each m = 2..{v - 1}.\n")
    f.write("# Format: m: n_1 n_2 ... n_K   (each line re-verified with exact rational arithmetic)\n")
    for m in range(2, v):
        f.write(f"{m}: " + " ".join(map(str, wit[m])) + "\n")
sys.exit(0 if bad == 0 and not contains_v and (expected is None or total == expected) else 1)
