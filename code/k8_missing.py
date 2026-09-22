#!/usr/bin/env python3
"""
k8_missing.py -- read the merged k = 8 bitmap and list the integers missing at k = 8.
Jude Wallis, 2026-09-21.

  k8_missing.py BITMAP.gz A B            print every m in [A, B) that occurs in no 8-term representation
  k8_missing.py BITMAP.gz --check FILE [A B]
                                         compare FILE (one integer per line, '#' comments allowed) with
                                         the missing set on [A, B); default range is [2, max(FILE) + 1)

Bit m of the bitmap (byte m >> 3, bit m & 7) is set iff m occurs as a denominator in at least one
representation 1 = 1/n_1 + ... + 1/n_8 with n_1 < ... < n_8. It is the OR of the nine ef_enum
shard bitmaps and covers every m < 2^26.

By the nesting lemma, anything missing at k = 9 is missing at k = 8, so the missing set here is
the full candidate list for v(9) below 2^26.
"""
import gzip, sys

bm = gzip.open(sys.argv[1], "rb").read()
LIM = len(bm) * 8
seen = lambda m: (bm[m >> 3] >> (m & 7)) & 1

if sys.argv[2] == "--check":
    got = sorted({int(x) for x in open(sys.argv[3]) if x.strip() and not x.lstrip().startswith("#")})
    lo, hi = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (2, max(got) + 1)
    got = [m for m in got if lo <= m < hi]
    want = [m for m in range(lo, hi) if not seen(m)]
    extra, absent = sorted(set(got) - set(want)), sorted(set(want) - set(got))
    print(f"range [{lo}, {hi}): expected {len(want)} missing integers, file has {len(got)}")
    print(f"  in file but occur at k = 8: {len(extra)} {extra[:10]}")
    print(f"  missing at k = 8 but not in file: {len(absent)} {absent[:10]}")
    sys.exit(1 if extra or absent else 0)

A, B = int(sys.argv[2]), int(sys.argv[3])
if B > LIM:
    sys.exit(f"bitmap only covers m < {LIM}")
for m in range(max(A, 2), B):
    if not seen(m):
        print(m)
