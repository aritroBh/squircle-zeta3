import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
result = r'''# W5C2 — EXACT certificate for M4's weight-5 part: W5 = N5 + R5_min3 (Kimi K3)

## Verdict: EXACT / COMPUTER-ASSISTED — all 8 systems verified exactly.

For D = W5 − (N5 + R5_min3) as exact level-8 H-word vectors (D_rat, D_sq = rational and √2
components, weight 5, 46656 admissible words), we produced, for each of the four Galois
sectors s ∈ {(+,+),(+,−),(−,+),(−,−)} and each component X ∈ {rat, sq}, an explicit
rational vector y^(s,X) with

  Σ_i y^(s,X)_i · R^(s)_i  =  Φ_s(D_X)   EXACTLY (Fraction arithmetic, residual 0 on every sector column),

where R^(s)_i are the W5C1 sector rows (proved relation folds). Together with the exactly
verified fold/unfold lemma and the exact unfolded assembly, this proves D = 0 as a value,
i.e. **W5 = N5 + R5_min3** (given the rows are proved relations — the planner's premise).

## The exact statements verified (all by my own code, printed in logs)

1. **Trimming (mod p1 = 16777213)** — LeanEchelon sparse GE per sector; S = independent rows,
   C = pivot columns, D_rat and D_sq reduce to 0:
   | sector | rows | cols | rank = \|S\| = \|C\| | S families | quotient |
   |---|---|---|---|---|---|
   | (+,+) | 83147 | 12218 | 12156 | engine 12120 + framew5 36 | 62 |
   | (+,−) | 76044 | 11164 | 11101 | engine 11065 + framew5 36 | 63 |
   | (−,+) | 76181 | 11164 | 11041 | engine 11041 (all) | 123* |
   | (−,−) | 82557 | 12110 | 11999 | engine 11999 (all) | 111* |
   (* quotient to the found span; D was already in span at that rank — full rank not needed.)
   Rank matches W5C1's Wiedemann lower bounds where rank is full (pp 12156 = minpoly 12156).
2. **Second prime (p2 = 16777199)**: rank(S) = |S|, rank(R_{S,C}) = |S| (det ≠ 0), D ∈ span(S) — all 4 sectors PASSED.
3. **Exact solve** y·R_{S,C} = D_C over Q (custom Dixon lifting mod p1, both RHS at once, LU
   via flint nmod_mat.inv, early-terminating rational reconstruction, stability repeat):
   | sector | Dixon steps | y_rat height (num/den digits) | y_sq height | A·Y = B mod p2 |
   |---|---|---|---|---|
   | (+,+) | 577 | 2078 / 2065 | 2081 / 2069 | ✓ |
   | (+,−) | 563 | 2023 / 2012 | 2027 / 2018 | ✓ |
   | (−,+) | 553 | 1991 / 1981 | 1991 / 1982 | ✓ |
   | (−,−) | 552 | 1990 / 1979 | 1989 / 1980 | ✓ |
4. **Exact verification (step_04, Fractions)**: Σ y_i R_i = Φ_s(D) on ALL sector columns —
   **residual 0 for all 8 systems**. Preimage check: every used sector row = (1/λ)·fold of its
   unfolded source row (recomputed from rels_*.pkl / validated_w5.pkl), λ matches.
5. **Fold/unfold lemma**: u = Σ_s (1/4) Σ_{g∈G} χ_s(g) g(Φ̃_s(u)) verified on ALL 46656 unit
   vectors (exact, all 4 sectors). Σ_s n_s = 46656 = N ✓.
6. **Exact unfolded identity (step_05, Fractions)**: D = Σ_s Σ_i (y_i λ_i/4) Σ_g χ_s(g) (g·r_i)
   holds on ALL 46656 words, **residual 0**, for both D_rat and D_sq.
7. **Standalone verify_W5.py** (fresh process, reads only certificates + W5C1 source files):
   orbit structure OK; lemma OK; all 8 sector equations residual 0; unfolded assembly over all
   4 sectors mod 2^61−1 residual 0 (rat 25406 words, sq 28535). **VERDICT: EXACT.**

## Why this proves W5 = N5 + R5_min3

The relation rows are vectors r with r·V = 0, where V is the vector of level-8 weight-5
H-values (proved relations: see row families below). The exact unfolded identity gives
D = Σ c_{g,r} (g·r). For each term, (g·r)·V = Σ_w r[w] V[gw] = r·(gV) = g(r·V) = g(0) = 0,
because (i) r·V = 0 (proved relation), (ii) g ∈ Gal(Q(ζ8)/Q) = {1, τ: j→−j, σ: j→3j, στ}
acts on the values by conjugation of the cyclotomic periods, and (iii) r has rational
coefficients. Hence D·V = 0, i.e. W5 = N5 + R5_min3 as real numbers.
No Galois-closure of the row SET is needed (step_06: 36016/40392 transforms are literally in
the files anyway; the rest are covered by the conjugation argument above).

## Row families used by the certificate (for the planner's provenance check)

- **engine/standard (std)**: regularised double shuffle (Ihara–Kaneko–Zagier) + regularised
  distribution d ∈ {2,4,8} (Zhao) — used in all 4 sectors (pp 9835, pm 9019, mp 8980, mm 9681 rows).
- **engine lifts**: shuffle lifts of proved lower-weight relations — lift1 (2248/2007/2037/2275),
  lift2old (20/8/12/22), lift3old (4/6/6/8), lift4old (6/11/5/10).
- **frame-w5** (288 weight-5 frame double-shuffle rows in real Möbius frames, validated by W5C1
  to 2.8e-113): **36 rows used in (+,+), 36 in (+,−), 0 in (−,+) and (−,−)** (engine alone
  spans D in the (−,·) sectors). This is consistent with W5C1's certified necessity of
  frame-w5 in (+,+) and (+,−).
- **frame-234 lifts**: NOT used in any sector.

Used rows are identified in the certificates by (family, index-in-family-file, label), so the
planner can regenerate each from rels_std/lift1/lift2old/lift3old/lift4old.pkl and
validated_w5.pkl and confirm they are proved relations.

## Files

- `certificate_{pp,pm,mp,mm}.json` (~90 MB each): sector, rank, |S|, primes, Dixon stats,
  heights, family counts, and per used row {S_index, fam, famidx, tuple, y_rat, y_sq,
  src {ufam, uidx, label, λ}}. Σ_i y·tuple = Φ_s(D) exactly (both components).
- `verify_W5.py` — single standalone verifier (rerun: full pass, ~60 s, exit 0).
- `work/`: step01/02_*.pkl (S, C, checks), yexact_*.pkl (exact y), feed_*.pkl, foldmap_*.pkl.
- Steps: step_01b (trim), step_02 (p2 checks), step_03 (Dixon solve), step_04 (exact verify +
  certificates), step_04a (fold replay/preimages), step_05 (exact unfolded identity),
  step_06 (Galois set-closure info), w5c2_lib/fold/dixon (modules), progress.log, summary.md.

## Timings (printed by the code)

| step | pp | pm | mp | mm |
|---|---|---|---|---|
| step_01b trim | 679 s | 992 s | 385 s | 371 s |
| step_02 p2 checks | 249 s | 251 s | 232 s | 229 s |
| step_03 exact solve | 963 s | 805 s | 776 s | 905 s |
| step_04 verify+cert | 20 s | 19 s | 16 s | 16 s |
+ step_04a 34 s, verify_W5 60 s, step_05 192 s. nmod inv(12156) = 406 s (Dixon LU).
Total pipeline compute ≈ 2.3 h (plus development/benchmarking).

## Honest caveats

- The certificate depends on the correctness of W5C1's D (D_rat.json/D_sq.json — validated
  numerically to 1e-110 there) and on the row families being proved relations (planner's
  premise; engine families are the cited theorems; frame-w5 rows are W5C1's validated set).
- The Galois conjugation of the cyclotomic values (step "Why this proves") is the standard
  action on convergent/regularised iterated integrals at roots of unity; it is the same
  action W5C1's sector construction is built on.
'''
open('result.md','w').write(result)
summary = r'''# W5C2 — exact certificate for W5 = N5 + R5_min3 (status: COMPLETE — EXACT, all 8 systems)

DONE. All 4 sectors × {rat,sq} verified EXACTLY (Fractions, residual 0) + exact unfolded
identity on all 46656 words + standalone verify_W5.py full pass. Heights ~1990-2081 digits.
Families used: engine std + lift1..4 + framew5 (36 in pp/pm, 0 in mp/mm); NO frame234.
See result.md. Verifier: python verify_W5.py (exit 0, ~60 s).
'''
open('summary.md','w').write(summary)
with open('progress.log','a') as f:
    f.write('2026-10-29 06:40:00 result.md + summary.md written. W5C2 COMPLETE: EXACT certificate for W5 = N5 + R5_min3, all 8 systems verified exactly.\n')
print('result.md and summary.md written')
print(subprocess.run(['bash','-c','ls -la'], capture_output=True, text=True).stdout[:2500])