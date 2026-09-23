# erdos293-check

An independent check of v(k) for [Erdős Problem #293](https://www.erdosproblems.com/293), together with a
repaired proof of the nesting lemma that the published lower bound on v(9) relies on.

v(k) is the least integer greater than 1 that never occurs as a denominator in any representation
1 = 1/n₁ + … + 1/n_k with n₁ < … < n_k. Erdős and Graham asked how fast it grows.

By Jude Wallis, 2026-09-19. Builds on computations by papanokechi, jonathanwxh-cell and
Muhamadiev Faridun, and on van Doorn–Tang, [arXiv:2512.22083](https://arxiv.org/abs/2512.22083).
Methods note: I used an AI coding assistant while writing the programs in `code/`; I ran them and
verified the arguments myself.

## Results

| | |
|---|---|
| v(3..8) | **4, 11, 17, 103, 733, 27539**, by full enumeration of every representation |
| representations found, k = 3..8 | 1, 6, 72, 2320, 245765, **151,182,379** — equal to [A006585](https://oeis.org/A006585), so nothing was missed |
| positive half of v(8) | an explicit representation containing m for **every** 2 ≤ m < 27539 (`data/ef.k8.witnesses.txt`, 27,537 of them), each re-verified with exact rational arithmetic |
| negative half of v(8) | 27539 occurs in no 8-term representation, confirmed twice: by the full enumeration, and by an independent per-candidate search using a different algorithm (`code/percandidate.py`) |
| integers missing at k = 8 | the complete list below 200,000: 5,466 of them, beginning 27539, 32411, 33647, 34919, 35279, 37199, … (`data/k8_absent_integers_below_200000.txt`). 1,572 are composite, the first being 53294 |
| all of k = 8 | every denominator below 2²⁶ that occurs at k = 8, as a bitmap (`data/k8_present_below_2pow26.bitmap.gz`). 57,146,348 integers in that range are missing, and by the nesting lemma only these can be missing at k = 9. The 1,157,643 missing below 3,000,000 are also written out (`data/k8_absent_integers_below_3000000.txt.gz`) |
| A097048 cross-check | Hugo van der Sanden's table of the least number of distinct unit fractions for (p−1)/p, p prime ≤ 800399 (linked from [A097048](https://oeis.org/A097048)), agrees with the enumeration for every odd prime and every k = 3..8: p occurs in a k-term representation of 1 exactly when his count is at most k−1 (`logs/hvds-crosscheck.log`) |
| v(9) | **13,856,993**, see below |
| nesting lemma | see below |

## v(9) = 13,856,993

Same two halves as v(8).

Every m with 2 ≤ m ≤ 13,856,992 has an explicit certificate: nine distinct denominators, one of them m,
whose reciprocals sum to exactly 1. These don't depend on any search being complete.

13,856,993 is in no 9-term representation. `code/findm_complete.c` goes through every way to write
1 − 1/13856993 as 8 distinct unit fractions, none of them 1/13856993. It ran as 2000 shards, found
nothing, and skipped no branch (`logs/v9-missing-13856993.log`, 4.07·10¹¹ nodes, about 30 min on a
10-core laptop). Before that it was run on the known cases: 474 representations with 8 terms contain
27538 (`percandidate.py` also gets 474), and none contain 27539. At k = 7 and k = 6 it finds none for
733 and 103. The full table is at the top of the log.

```bash
clang -O2 -o findm_complete code/findm_complete.c -I$(brew --prefix gmp)/include -L$(brew --prefix gmp)/lib -lgmp
./findm_complete 27538 3 1 0 8 | tail -1          # 474 solutions, skips 0
for w in $(seq 0 1999); do ./findm_complete 13856993 6 2000 $w 9 | tail -1; done   # parallelize this
```

The certificates are attached to the
[v9-certificates release](https://github.com/EconLearn/erdos293-check/releases/tag/v9-certificates):
`v9_certificates_2_to_3000000.txt.xz` (65 MB), `v9_certificates_3000001_to_10000000.txt.xz` (113 MB) and
`v9_certificates_10000001_to_20000000.txt.xz`. The last file has gaps above 13,856,992, listed at its
end. Some certificates in it are ntysdd's, from the solver they posted on #403. To check:

```bash
xz -dk v9_certificates_*.txt.xz
python3 code/check_v9_certs.py v9_certificates_2_to_3000000.txt v9_certificates_3000001_to_10000000.txt v9_certificates_10000001_to_20000000.txt 2 13856992   # ~10 min
```

It re-adds every line with `fractions.Fraction` and confirms that no m in the range is missing.

How the certificates were found (counts in `logs/v9-certificates.log` for m ≤ 3,000,000 and
`logs/v9-certificates-3M-10M.log` above that):

- `code/findm.c`: depth-first search for (m−1)/m as 8 distinct unit fractions, none equal to 1/m,
  with a cap on the work per m. This settled 2,988,471 values.
- `code/split_stage.py`: if s divides m² and d = m − s occurs in an 8-term representation W, then
  1/d = 1/m + 1/(dm/s), so W − {d} + {m, dm/s} works when dm/s is not already in W (7,576 values).
- `code/lift.py`: an 8-term representation containing m, lifted by the nesting construction (1,155).
- `code/threesplit.c`: split one term t < m of a 7-term representation into three,
  1/t = 1/m + 1/x + 1/y, using all 245,765 seven-term representations
  (`data/k7_all_representations.txt.gz`). This settled all 2,797 values the others missed below
  3,000,000, in under a second, and 6,999,957 of the 7,000,000 integers from 3,000,001 to 10,000,000 in
  about a minute.
- `code/padic_jobs.py` + `code/findfrac.c`, for primes p with p − 1 missing at k = 8. The other multiples
  of p in the representation have to cancel the 1/p p-adically. The two simplest ways: one more multiple
  p(jp − 1), with the remaining 7 terms summing to 1 − j/(jp − 1), or two, pk₁ and pk₂ with
  (k₁ + 1)(k₂ + 1) = 1 + ℓp, with the remaining 6 terms summing to 1 − ℓ/(k₁k₂). Each (j) or (ℓ, k₁, k₂) is a smaller search
  (all 25 values below 10,000,000 that the other methods missed, the first being 6,156,713). Running
  each subproblem to completion rather than under a small work cap mattered: ntysdd's certificates for
  9121003, 9243467, 9305819 and 9307751 on #403 showed that the one-multiple family reaches primes the
  capped search had given up on.

Muhamadiev Faridun (forum, 14 Sep: v(9) > 2,108,538) and ntysdd
([#403](https://github.com/teorth/erdosproblems/issues/403): v(9) > 3,000,000 by a scan, with
certificates for 43,538 values) had the earlier bounds. These certificates were found separately; of
ntysdd's 43,538, only 118 coincide with the ones here.

## The nesting lemma

If m occurs in a k-term representation (k ≥ 3), it occurs in a (k+1)-term one. This is what makes the
sets of missing integers nested, hence v(k) non-decreasing, and it is what the bound v(9) ≥ 27539 uses.

The usual one-line argument — split the largest denominator N into N+1 and N(N+1) — **does not cover the
case where m is itself the largest denominator**. A repair, with the case analysis in `MATH-NOTES.md` §4:
split the *second* largest denominator a as 1/a = 1/(a+d) + 1/(a + a²/d) for a divisor d of a² with d < a.
Both new denominators exceed a, so only m can collide; the pairs for different d are disjoint, so m blocks
at most one d; a composite a has at least two choices; and a prime a = p is impossible by looking at the
power of p in the sum.

`code/nesting_lemma_check.py` applies exactly this construction to every (representation, m) pair with
k ≤ 7 — 1,720,355 pairs at k = 7, of which 245,765 are the awkward case — and verifies each result exactly.

## Reproducing

```bash
clang -O2 -o code/ef_enum code/ef_enum.c
mkdir -p out && ./code/ef_enum 7 0 1 out/ef            # a few seconds
python3 code/merge_and_verify.py out/ef 7 1 245765     # count, v(7), and every witness re-verified

python3 code/percandidate.py vk 7                      # the second method, independently
python3 code/percandidate.py count 8 27539             # ~30 s: exhaustive, finds nothing
python3 code/nesting_lemma_check.py 6                  # the lemma, constructively

# candidates for v(9): integers in [A, B) missing at k = 8; or test a search program's k = 8 output
python3 code/k8_missing.py data/k8_present_below_2pow26.bitmap.gz 2108539 3000000
python3 code/k8_missing.py data/k8_present_below_2pow26.bitmap.gz --check your_k8_missing.txt 2 200000
```

k = 8 is the same program run as nine shards (`./code/ef_enum 8 <w> 9 out/ef` for w = 0..8), about
90 minutes wall-clock on a 10-core laptop, then `merge_and_verify.py out/ef 8 9 151182379`. Console
output of the runs behind every number above is in `logs/`.

## Files

- `code/ef_enum.c` — full enumeration. Exact integer arithmetic; the last two denominators are solved
  through the divisors of q² rather than by search.
- `code/percandidate.py` — the independent method: fix m, subtract 1/m, solve the rest; the two-term case
  uses a coprime-pair parametrization instead. Shares no code with the above.
- `code/merge_and_verify.py` — merges shard output, recomputes v(k), and re-adds every witness with
  `fractions.Fraction`.
- `code/nesting_lemma_check.py`, `code/check_nesting.py` — the lemma, constructively and from the bitmaps.
- `code/k8_missing.py` — reads the k = 8 bitmap: lists the missing integers in a range, or compares a
  list produced by another program with them.
- `code/compare_hvds.py` — the A097048 cross-check above.
- `code/findm.c`, `code/split_stage.py`, `code/lift.py`, `code/threesplit.c`, `code/padic_jobs.py`,
  `code/findfrac.c`, `code/padic_build.py` — the v(9) certificate search; `code/check_v9_certs.py` — the
  standalone checker for the certificate files.
- `code/findm_complete.c` — exhaustive search for one m, sharded; how 13,856,993 was ruled out.
- `data/k7_all_representations.txt.gz` — every 7-term representation of 1 (245,765 = A006585(7)).
- `data/k8_present_below_2pow26.bitmap.gz` — bit m (byte m >> 3, bit m & 7) is set iff m occurs in some
  8-term representation; the OR of the nine shard bitmaps. sha256 of the uncompressed 8 MiB file:
  `264d491519f0a671cc9cdb394370e86b3586972c95cafe4405e2fcdd48e9946f`.
- `MATH-NOTES.md` — the method and the proofs, written out.
