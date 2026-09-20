"""check_nesting.py -- is every denominator that occurs with exactly k terms also one that occurs
with exactly k+1 terms?  (Needed for 'v(k) is non-decreasing' and for the 'at most k terms'
reading of the definition.)  Uses the bitmaps written by ef_enum (denominators < 2^26 only).
Jude Wallis, 2026-09-19. Written with the help of an AI coding assistant; checked and run by me."""
import sys
def load(prefix, k, nw):
    bm = None
    for w in range(nw):
        b = open(f"{prefix}.k{k}.w{w}.bitmap", "rb").read()
        bm = int.from_bytes(b, "little") if bm is None else bm | int.from_bytes(b, "little")
    return bm
runs = {3: ("data/run/ef", 1), 4: ("data/run/ef", 1), 5: ("data/run/ef", 1), 6: ("data/run/ef", 1), 7: ("data/run/ef", 1)}
if len(sys.argv) > 1: runs[8] = ("data/run8/ef", 9)
bms = {k: load(p, k, nw) for k, (p, nw) in runs.items()}
for k in sorted(bms):
    if k + 1 in bms:
        missing = bms[k] & ~bms[k + 1]
        cnt = bin(missing).count("1")
        first = []
        m = missing
        while m and len(first) < 10:
            low = m & -m; first.append(low.bit_length() - 1); m ^= low
        print(f"k={k}: {bin(bms[k]).count('1')} denominators < 2^26 occur; of these {cnt} do NOT occur with k+1 terms", first)
