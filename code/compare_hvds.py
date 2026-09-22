#!/usr/bin/env python3
"""
compare_hvds.py -- compare Hugo van der Sanden's table of the least number of distinct unit
fractions needed for (p-1)/p (p prime, p <= 800399; 'results-single' in github.com/hvds/seq,
least_eg/, linked from OEIS A097048) with the ef_enum bitmaps for k = 3..8.
Jude Wallis, 2026-09-21.

  compare_hvds.py results-single PREFIX_K3_TO_K7 K8_BITMAP.gz

PREFIX_K3_TO_K7 is the output prefix of single-worker runs ./code/ef_enum k 0 1 PREFIX, k = 3..7.

m occurs in a k-term representation of 1 exactly when (m-1)/m is a sum of k-1 distinct unit
fractions none of which is 1/m. Hugo's count does not exclude 1/p, so the two can in principle
differ; the check is that for every odd prime p in the table and every k,
    [p occurs in some k-term representation of 1]  ==  [his count for p is at most k-1].
"""
import gzip, re, sys

hv = {}
for line in open(sys.argv[1]):
    m = re.match(r"(\d+): (\d+) ", line)
    if m:
        hv[int(m.group(1))] = int(m.group(2))
print(f"{len(hv)} primes in the table, largest {max(hv)}")

bad_total = 0
for k in range(3, 9):
    bm = gzip.open(sys.argv[3], "rb").read() if k == 8 else open(f"{sys.argv[2]}.k{k}.w0.bitmap", "rb").read()
    seen = lambda m: (bm[m >> 3] >> (m & 7)) & 1
    odd = [p for p in hv if p > 2]
    bad = [p for p in odd if bool(seen(p)) != (hv[p] <= k - 1)]
    bad_total += len(bad)
    print(f"k = {k}: {len(odd)} odd primes compared, {len(bad)} mismatches {bad[:10]}")
print("OK" if bad_total == 0 else f"{bad_total} MISMATCHES")
