(no FINAL answer; see transcript.json)

Last reply:
FINAL:

# TOW3R — referee report on TOW3 (tower law extended to depth 6)

**Verdict: (a) tower law at depths δ ≤ 6 — PROOF COMPLETE (MINOR wording fixes). (b) weight-5 law — PROOF COMPLETE in its logic chain; the two finite kernel checks (k=5 kernels 40/40; ±i simple-pole zeros at k=6..9) were NOT independently rerun by me (recorded in TOW3's transcript; classified "not run by referee", not a gap in TOW3). (c) Task D — unfinished, accurately stated.** Every load-bearing symbolic certificate was re-executed by me and prints remainder/numerator 0. The guess-and-verify for the fitted e(5,1), e(6,1), e(6,2) **closes**, by the same argument as at depth 4. One stale script (`step_33.py`) is buggy as committed — the corrected run exists only in the transcript; flagged, not load-bearing.

## 1. What I ran (all exact arithmetic; residual = 0 means *exactly* zero)

| # | check (my rerun / my own code) | result |
|---|---|---|
| 1 | `prove_V56.py` all 7 non-adjacent identities (5,0)..(6,3) | numerator **0** each (0.0–5.0 s) |
| 2 | `prove_adj56.py` all 7 adjacent identities | remainder·DEN **0** each: (5,0) 169 s, (5,1) 25 s, (5,2) 1.7 s, (6,0) 570 s, (6,1) 103 s, (6,2) 14 s, (6,3) 0.8 s |
| 3 | `prove_E3b56.py` δ=5,6 | numerator **0** (7.5 s, 15.3 s) |
| 4 | TOW2 base from TOW3's byte-identical copies: `prove_S12`, `prove_S34 3/4`, `prove_tops`, `prove_adj34` (3,0),(3,1),(4,0),(4,1),(4,2); TOW2 `step_48` V(3), `step_52` slope(4,1), `step_49` V(4) | all remainders **0** (V(4): 256 s) |
| 5 | `L_parsed[r] ≡ tow2lib.L(r)` **symbolically** in Q(i,√2)(d,m), r=1..6 (TOW3 had only 3 random points at r=5,6) | **equal** (my step_04; √2-aware reduction) |
| 6 | ρ(δ)/ρ(0) = a_δ vs TOWN1 towers, δ=3,4: 723 cells; δ=5,6: 720 cells each | **0 bad**; a₅,a₆ displayed closed forms ≡ Bell assembly symbolically |
| 7 | from-scratch A(v) log-derivatives (genintegrand.S_poly, no tow2lib): a₅,a₆ at 11 keys incl. (9,0), diagonal (26,26) d=52 | all match data ratios (residual 0) |
| 8 | actual off-diagonal Newton coefficients extracted from raw data vs closed forms: depth 5: **1932**, depth 6: **2476** coefficient-cells | **0 bad** (covers a_δ, e(5,1), e(6,1), e(6,2), zigzag tops) |
| 9 | actual diagonal ê vs ê-formula: 66 + 88 coefficient-cells | **0 bad** |
| 10 | empirical E3 from raw data: Δ_pp built from actual towers, fed via direct+Fold (TOW1's formula) | **= 0** for δ=5,6, p=3..10 (16/16) |
| 11 | **my own exact T1 recursion** (τ₁,τ₂ computed from the definition in exact Q(√2,i) arithmetic, landing in Q(i), = tow2lib table n≤7) vs data: degrees 12..30, k=d+1..52, depths ≤6 | **44443 fresh cells, 0 bad** |
| 12 | `check_tau(7)`; conj(τ1_n)=(−1)ⁿτ2_n and vice versa, n=0..7 | **True** |
| 13 | fit pickles ≡ tow3lib: e(4,1), e(5,1), e(6,2) (`e_fits.pkl`), e(6,1) (`e61_closed.pkl`) | **equal** (after `cancel`) |
| 14 | e00-ratio constants in the scripts vs data p=3,4,5: rB0d = 2i(4p−1)/(p+1), rB0 = (4p−1)(4p−3)/((p+1)p), rA0 = i(4p−3)/(p+1); e00_pp anti-real, e00_{p,p−1} real | **all confirmed** |

Digit counts: all checks are exact (Fraction/sympy Q(i,√2) arithmetic); residuals are exactly 0, no floating point anywhere.

## 2. Findings by item

**1. Genuinely symbolic — OK.** `prove_V56/adj56/E3b56` treat (d,m) resp. p as free symbols; identities are rational identities over Q(i) (remainder 0 ⇒ valid wherever denominators ≠ 0; all denominators are products of 2d−2j−1, 4p−odd, p, p+1 — never zero at integers in range). Fitted objects: e(5,1), e(6,2) (TOW2 fits), e(6,1) (TOW3 fit, 60 keys). **Closure:** TOW1's reduction theorem needs only *existence* of an e-family satisfying E0 + all N(δ′,s) (δ′≤6, both key types) + E3; the identity set at δ≤6 is complete (δ′≤4 inherited from TOW2 — byte-identical libs/scripts, all rerun green; δ′=5: s=0,1,2; δ′=6: s=0..3; adjacent same seven; E3 δ=5,6), and each fitted form enters several proved identities (e.g. e(6,1) enters (6,0),(6,1), adjacent (6,0),(6,1), E3(6)). Finite descending induction in degree at fixed k (initial data ρ at degree k, uniqueness of T1's solution) gives U=u. The fit provenance is inessential — exactly as at depth 4. `e_fits.pkl` also contains (7,·),(8,·) fits (incl. a degenerate (7,3) with b=0): **unused** at δ≤6 (Task D only). One GAP closed by referee: TOW3 verified the simplified L-cache against the proved L(r) only at 3 random points (r=5,6); my step_04 proves equality symbolically for r=1..6, so ρ(5),ρ(6) now rest entirely on proved ingredients.

**2. ρ(5), ρ(6) — OK.** I re-derived L_r = (log A)^⁽ʳ⁾(i)/r! by hand (A(v) = −(√2/2)8^d v^{2d−1}S_m(t(v))(v+i)^{2−2d}/((v²−a²)(v²−b²)), t(v)=v/4+1/2−1/(4v), t′(i)=0; Faa di Bruno blocks of size ≥2 with t⁽ʳ⁾(i)=(−1)^{r−1}r!/(4i^{r+1}); sl↔R conversion standard) — it is exactly tow2lib.L(r); Bell assembly a_δ = Σ_partitions ∏L_sz^mult/mult! = A⁽ᵟ⁾(i)/(δ!A(i)) = ρ(δ)/ρ(0). Verified: displayed N₅,N₆ forms ≡ Bell (symbolic); 720 cells each δ=5,6, 0 bad; 11-key from-scratch check. **MINOR FIX:** committed `step_33.py` is the *stale* version with `expected = (1+[p>q])·a_δ`; transcript shows its first run "tested 720 bad 696" and the corrected rerun (0 bad) was never written back. The correct normalization (confirmed by my checks): a_δ = ρ(δ)/ρ(0), the 1+[p>q] being inside ρ(0)=e00 — the writeup states this correctly; the script must match.

**3. Newton identities, fold, ê, E3 — OK.** All rerun green (table). I re-derived the adjacent feed: rB0 = rA0(2p,p)·rB0(2p−1,p) = e00_{p+1,p}/e00_{p,p−1} (real); rB0d = e00_{p+1,p}/e00_pp (with the (1+[p>q]) doubling); the alternating-sum conversion Σ_t (−1)^t ê(δ′,s+t); the fold exponents (−i)^{4p}i^{4p}=1; and — the CX8-recurring subtlety — **e00_pp is anti-real** (conj(e00)=−e00) so the actual conjugate in e00-units carries an extra minus, which is exactly `Dc = −Δ̄` in `prove_E3b56.py`; the adjacent child's e00 is real, no flip. The code's signs are correct; my empirical E3 from raw data (item 10) independently reproduces direct+Fold = 0 at δ=5,6, p=3..10, and the diagonal ê extraction (item 9) confirms the ê-map at 154 fresh coefficients.

**4. Logic chain to (b) — OK as stated.** WTS Lemma 3 (PROVED): T for degree n=d+1 at depths [0, 2d−2] for k≥n+1 ⇒ degree-d kernel residues polynomial in k of degree ≤d−1 for k≥d+2. For d=4 this needs degree-5 towers to depth 6 = 2·5−4 — exactly TOW3's theorem; degree 5 is odd ⇒ no diagonal key ⇒ law holds from k=5 ⊇ the k≥6 hypothesis. Lemma 2 (fixed R-linear tower→kernel map, PROVED) ⇒ W₅(k) cubic for k≥6. k=5 (=d+1) is outside Lemma 3; TOW3 settles it by the exact kern_5-vs-cubic(6..9) check (40/40) + cubicity k=6..13, and the ±i simple-pole polynomial (degree ≤3) vanishes for k≥6 given zeros at k=6..9. **Not run by referee:** the 40/40 kernel comparison and the ipos/ineg zeros (WTS kern_k.pkl) — TOW3's transcript records them green; repairable in minutes (load kern_5..kern_9, compare). Weights polynomial in k for ALL k≥w: w=2 (TOW1), w=3 (WT3+TOW1), w=4 (TOW2 + k=4 exact, per TOW2R/CX8), **w=5 (here, k=5 exact)**; w≥6 open (needs δ≥8). LEAD1/LEAD2 top table: zigzag-top *equations* proved all s (non-adjacent; adjacent tops only s≤3, now incl. (5,2),(6,3)); numerators A_{s+1} are data-verified (TOW2 s≤6; 104 keys here) and become the true coefficients a posteriori where the tower law holds (δ≤6) — TOW3 states this accurately.

**5. Independent data — OK.** My own exact T1 (item 11): 44443 fresh cells (degrees 12..30, k≤52, δ≤6), 0 bad — the check TOW2R explicitly left unrun. Plus items 6–10 (~7500 coefficient/rho cells), all exact.

**6. Task D — noted only.** Claims: Res_{−a} map derived and validated vs SUB1 (d≤5; kernel spot checks); M_d reduced to one sub-top value per d (validated d=2..7 = SUB1's −1,10,15,131,394,3151); Catalan pattern for e(2s,s−1) CONJECTURE; M_d closed form NOT found. Accurately stated as unfinished; not part of the verdict.

## 3. Wording fixes (MINOR)

1. Restore the corrected `step_33.py` (drop the `mult` factor; transcript run 2). 2. proof_depth6.md §1: replace "exact random-point check at r=5,6" by symbolic equality (supplied: my step_04, r=1..6). 3. §3: state the anti-real e00_pp conjugation sign (`Dc=−Δ̄`) and the reality of e00_{p,p−1} (CX8 #2 recurs). 4. State domains: depth δ needs 2d−4≥δ (δ=5,6 ⇒ d≥5 — not d≥6); adjacent (p,p−1) needs 4p−6≥δ (p≥3); k≥d off-diagonal, k≥2p+1 diagonal stable, k=2p exceptional. 5. Mark `e_fits.pkl` entries (7,·),(8,·) as unused candidates (Task D); the (7,3) fit is degenerate. 6. TOW3's `step_48.py` writes `prove_adj56.py` (code-gen; crashes read-only) and `step_52.py` is a driver — keep TOW3 step numbers out of proof citations (currently OK). 7. Display the e00-ratio constants rA0/rB0/rB0d in the text (verified vs data here). 8. `verify_tow3.py` still absent — say so (per-lemma certificates suffice; CX8 #8 recurs).

## 4. Deliverables in `new_results/TOW3R_review/`

`progress.log`; `step_02.py`–`step_17.py` (+`.out`) for items above; `tow3_copy/` (byte-verified copies). Key evidence: step_04 (L symbolic, r=1..6: all True); step_05 (ρ: 723/723, 720/720 ×2, 0 bad; Bell≡displayed); step_06/07 (11 keys from-scratch A(v): all True); step_08/09a–d/10/11 (all prove_* reruns green); step_12 (fits≡lib; 1932+2476 cells 0 bad); step_13 (ê 66+88; E3 empirical 16/16 zero); step_14 (T1: 44443 cells, 0 bad); step_15 (τ rules); step_16–18 (δ≤4 base green).

`verify_tow3r.py` (< 10 min subset; slow certificates listed with exact command lines):
```python
# verify_tow3r.py -- run in TOW3R_review/ with project venv python. PASS/FAIL per item.
import subprocess, sys, os, time
W = os.path.dirname(os.path.abspath(__file__)); PY = sys.executable
def run(name, script, *args, cwd=None, timeout=1700):
    t0=time.time()
    r = subprocess.run(['nice','-n','10',PY,script,*map(str,args)], capture_output=True, text=True,
                       cwd=cwd or W+'/tow3_copy', timeout=timeout)
    ok = r.returncode==0
    print(f'{"PASS" if ok else "FAIL"}  {name}  ({time.time()-t0:.0f}s)')
    return ok
ok = True
ok &= run('L_parsed==L(r) symbolic r=1..6', W+'/step_04.py')          # ~10 s
ok &= run('rho(3..6) 2886 cells + Bell==displayed', W+'/step_05.py')  # ~20 s
ok &= run('tau table + conjugation rules', W+'/step_15.py')           # ~5 s
ok &= run('fits==lib; Newton coefs 1932+2476 cells', W+'/step_12.py') # ~15 s
ok &= run('diagonal ehat; empirical E3 (16 cases)', W+'/step_13.py')  # ~60 s
ok &= run('own T1 recursion 44443 fresh cells', W+'/step_14.py')      # ~30 s
for s in [(5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3)]:
    ok &= run(f'prove_V56 {s}', 'prove_V56.py', *s)                   # ~8 s total
ok &= run('prove_adj56 (5,2)', 'prove_adj56.py', 5, 2)                # fast representative
ok &= run('prove_adj56 (6,3)', 'prove_adj56.py', 6, 3)
ok &= run('prove_E3b56 5', 'prove_E3b56.py', 5)
ok &= run('prove_E3b56 6', 'prove_E3b56.py', 6)
ok &= run('prove_S12', 'prove_S12.py'); ok &= run('prove_S34 3','prove_S34.py',3)
ok &= run('prove_S34 4', 'prove_S34.py', 4); ok &= run('prove_tops','prove_tops.py')
# SLOW (documented, rerun green in this review; command lines):
#   prove_adj56.py 5 0 (169 s) | 5 1 (25 s) | 6 0 (570 s) | 6 1 (103 s) | 6 2 (14 s)
#   TOW2 base: step_48 V(3) 30 s, step_52 slope(4,1), step_49 V(4) 256 s, prove_adj34 all five
# NOT RUN BY REFEREE: Task C finite kernel checks (kern_5 vs cubic(6..9) 40/40; ipos/ineg k=6..9)
print('ALL PASS' if ok else 'SOME FAIL')
```

**Bottom line:** the tower law at δ≤6 (all keys, all d, all k; off-diagonal k≥d, diagonal k≥2p+1 with k=2p exceptional) and, given the recorded finite kernel checks, W₅(k) = W₅″ + (k−5)E1 + C(k−5,2)E2 + C(k−5,3)E3 for every k≥5, stand as proved. No mathematical errors found; the fixes above are presentational/reproducibility issues.