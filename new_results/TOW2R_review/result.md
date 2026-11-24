(no FINAL answer; see transcript.json)

Last reply:
FINAL:

# TOW2R — referee report on TOW2 (tower law extended to depth 4)

**Verdict: PROOF COMPLETE for the tower law at depths δ ≤ 4 (all keys, all d, all k), with MINOR wording fixes.** The guess-and-verify step for the fitted e(4,1) **closes** (analysis below). Every symbolic certificate in TOW2 was re-executed by me and prints remainder 0. Consequence (b) (weight-4 law for all k ≥ 4) follows as stated. Partial results (c) are accurately stated. One committed script (`step_47.py`) is broken and `tow2lib.e_desc` is mis-documented — neither is load-bearing. Items I could not finish are listed as **not run by referee** (not as gaps in TOW2).

## 1. What the proof is and why it closes

TOW1's reduction theorem (PROVED, TOW1R/CX7 errata recorded): tower law at depths ≤ δ ⟺ E0 + the triangular Newton identities N(δ′,s) (0 ≤ s ≤ ⌊δ′/2⌋, δ′ ≤ δ) for every key + blindness E3(δ′). The identities at general δ are derived in TOW1 §5 (Newton substitution of u = Σ_s binom(k−d,s)e(δ,s) into T1; diagonal feeds via the alternating ê-conversion; fold via conjugate g_μ). For δ ≤ 4 the required NEW certificates are exactly:

| identity | script | my rerun |
|---|---|---|
| S′/S, S″/S at t* (all (d,m)) | prove_S12.py | remainder **0**, 1.6 s |
| S‴/S, S⁗/S at t* (all (d,m)) | prove_S34.py 3 / 4 | numerator **0**, 2.4 s / 6.8 s |
| top eqs, all s, non-adjacent | prove_tops.py | all brackets **0**, 1.2 s |
| V(3) = N(3,0), non-adjacent | **step_48.py** (mis-cited in text) | remainder **0**, 39.3 s |
| V(4) = N(4,0), non-adjacent | **step_49.py** (mis-cited) | remainder **0**, 335.5 s |
| slope N(4,1), non-adjacent | step_52.py | remainder **0**, 1.7 s |
| adjacent (3,0),(3,1),(4,0),(4,1),(4,2) | prove_adj34.py δ s | all **True** (2.4, 1.7, 20.5, 6.9, 100 s) |
| blindness E3, δ = 3, 4 | prove_E3b34.py | both **True** (4.0, 7.8 s) |

Diagonal keys need no separate identities (parent (p+1,p) has degree 2p+1, so the DL1 feed lands directly in the binom(k−2p−1,s) basis; the k = 2p exception is exactly E3). The set above is complete, hence by TOW1's theorem + uniqueness of T1's level-by-level solution, U = u: the tower law holds at δ ≤ 4.

**The e(4,1) guess-and-verify closes.** The fitted form enters *every* identity of the triangular system at δ = 4 — N(4,0) (step_49), N(4,1) (step_52), the top N(4,2) (prove_tops, independent of the numerator), all adjacent identities (prove_adj34), and E3(4) (prove_E3b34) — and each is an exact rational identity in (d,m) resp. p over Q(i), remainder 0 (rerun above). The reduction theorem requires only *existence* of an e-family satisfying the identities; verification of the complete set therefore proves the law regardless of how the candidate was found. (The fit is even inessential: e(4,1) is determined by the proved top e(4,2) and the parents through the functional equation N(4,1); the fit just solves it.) The zigzag numerator at depth 4 (A₃ = 2) is then the true top coefficient since U = u.

## 2. Independent checks I ran (not from tow2lib)

- **S-derivatives from the definition** (genintegrand.S_poly, exact Q(i) arithmetic): S^{(l)}(t*)/S(t*) vs N_l/Π_l — l = 1,2: **86 cells** (d = 4..12, all keys incl. q = 0), **0 bad**; l = 3,4: **86 cells**, **0 bad**. This validates the closed forms *and* the sum-formula bridge of prove_S12/S34 at every boundary shape.
- **Stirling tails** P_r, r ≤ 4: recurrence P_r(x+1) = (x+1−r)P_{r−1}(x)+P_r(x) symbolic **True**; anchors **True**; P_r(n) = S(n,n−r) for n = 0..29 (genintegrand.stirling2) **all equal**; P_r divisible by x^underline(r) (so the summand is polynomial after cancellation — the q = 0 boundary is exact).
- **First-principles Φ_p**: I re-derived the integrand factor from the master square root (−1/2[ε^d]Φ convention): Φ_p = −(√2/2)8^d4^{−D}v^{2d−1−D}S̃(v²+1)^{2−2d}/((v²−a²)(v²−b²)) — matches buildintegrand.py and TOW2's logA term by term (v-exponent 2d−1−D = −1, (2−2d)log(v+i), D = deg S = 2d confirmed for d = 2..12). My own Laurent extraction at v = i (QSI series inversion) reproduces binteg initial data for key (4,0) at δ = 0..4 (**5/5 exact**, e.g. ρ(0) = 5√2i/16).
- **Top-equation algebra**: I redid the even/odd top reductions by hand: E0 bracket and τ1₁-feed bracket (M¹ = (3i/2)(2d−1)) correct; the telescoping Π_s(d+1)/Π_s(d) = (2d−1)/(2d−2s−1) is an exact cancellation (proved in the writeup, correct); the odd-top identity −(X−2)/(2X−1) = −(X−1)/(2X−2) + (3/2)/(2X−1), X = d−s, is exact.
- **Fold exponents** in prove_adj34/prove_E3b34: (−i)^{4p−4−t} = (−i)^{−4−t}, i^{4p−6−δ} = i^{−6−δ} for integer p — legitimate; the session's earlier "False" remainders were (i) a structural `rem == 0` check and (ii) unreduced symbolic i-powers, not mathematical failures (transcript forensics; final scripts rerun green).

## 3. Findings by review item

1. **Genuinely symbolic — OK.** d, m (S-lemmas, V(3), V(4), slope), (d,s) (tops), p (adjacent, E3) are symbols; all summation indices are concrete with ranges depending only on the fixed order l or depth. The only fitted object, e(4,1), is guess-and-verify and closes (§1). R₃, R₄ were fitted then **proved** (prove_S34, rerun green) — §5 of proof_all_depths.md ("proof pending") is stale.
2. **ρ(3), ρ(4) — OK.** They depend only on sl₁, sl₂ (proved) via Faa di Bruno/Bell (t′(i) = 0 ⟹ only blocks ≥ 2). My independent S-checks (172 cells) + Φ re-derivation support the ingredients; TOW2's 72-cell ρ check and step_77's 1348-cell law check were **not rerun by me** (time) — flagged, not a gap in TOW2.
3. **Newton identities — OK** (table above; all rerun green). My own re-derivation of one adjacent identity was only partially completed (fold exponent and ê/rB0 conventions verified; full term-by-term reconstruction not finished) — honest note.
4. **Logic chain to (b) — OK as stated.** Weight w needs degree-w towers to depth 2w−4: w = 3 (TOW1+WT3, k ≥ 3), **w = 4 (here): degree-4 towers to depth 4 + WTS Lemma 3 ⟹ K(3) polynomial for k ≥ 5; k = 4 is exact data (WT4/WTS), not reproved — TOW2 says exactly this.** The ±i simple-pole coefficient vanishes identically once zero at 3 values k ≥ 5 (Lemma 3 criterion; exact data — not rechecked by me). W4(k) = W4″ + (k−4)D1 + C(k−4,2)D2 then holds for all k ≥ 4, upgrading paper thm:w4k (currently 4 ≤ k ≤ 52). w ≥ 5 correctly stated OPEN (needs δ ≥ 6). LEAD1: zigzag-top equations proved all s, numerator forced only s ≤ 2 — accurately stated.
5. **Independent fresh-data tower check (n = 12..20, k ≤ 52): NOT RUN by referee** (session cut short after the Φ fix). TOW2's own checks (72 cells ρ; 1348 cells k ≤ 14; planner reran prove_tops) stand unchallenged but un-extended by me.
6. **Deliverables:** this report is result.md; `verify_tow2r.py` is inlined below (it encodes exactly the checks I executed piecewise; each was run, under nice, during the review — total ≈ 10 min, dominated by V(4)).

## 4. Wording fixes (MINOR)

1. proof_all_depths.md §3: V(3)/V(4) non-adjacent are **step_48.py / step_49.py**, slope (4,1) is step_52.py (text says "steps 51,52,55"; step_51 prints source).
2. proof_all_depths.md §5: R₃, R₄ are **proved** by prove_S34.py (rerun: numerator 0 at l = 3,4) — delete "proof pending".
3. **tow2lib.e_desc is mis-documented and non-terminating**: the descent e(δ,s) = Λ(δ,s) − e(δ,s+1) recurses into e_desc_at(δ,s,d+1,·) (same (δ,s), degree +1) — RecursionError for every s < ⌊δ/2⌋ (confirmed by me for (2,0),(3,0),(4,0),(4,1)). Consequently **step_47.py as committed crashes** (the transcript shows the same crash). Not load-bearing (proofs use closed forms directly), but the result.md row "e(4,1) ... proved to satisfy all its identities" should cite step_49/step_52/prove_adj34/prove_E3b34, not the broken step_47. Remove or rewrite e_desc.
4. result.md: "ρ(δ) proved closed forms for δ ≤ 8" vs honest boundaries "ρ(5),ρ(6) assembly not finished" — reconcile (ingredients all proved; assembly/verification of ρ(5..8) pending; only ρ(3),ρ(4) are used at δ ≤ 4).
5. genintegrand.py docstring's v-space formula is stale (off by −(v²+1)/2 vs its own "coefficient × dt/dv" prescription, and a sign convention vs "−1/2[ε^k]Φ"); the code (buildintegrand.py) is correct — align the comment.
6. prove_tops.py: the telescoping step is comment-only; add the one-line induction check (Π_{s+1}(d) = Π_s(d)(2d−2s−3)) for completeness.
7. result.md: state explicitly that τ-constants were verified against the DL1/DL2 definitions (tow2lib.check_tau, n ≤ 7).

## 5. verify_tow2r.py (encodes the executed review; each piece ran green above)

```python
# verify_tow2r.py -- TOW2R consolidated referee check. Run in TOW2R_review/tow2_copy
# (copies of TOW2's scripts+libs) with project venv python. Prints PASS/FAIL per item.
import subprocess, sys, os, time
PY = sys.executable
def run(script, *args, timeout=1700):
    t0 = time.time()
    r = subprocess.run(['nice','-n','10',PY,script,*map(str,args)],
                       capture_output=True, text=True, timeout=timeout)
    ok = (r.returncode == 0)
    print(f'{"PASS" if ok else "FAIL"}  {script} {args}  ({time.time()-t0:.0f}s)')
    print('   ', r.stdout.strip().replace(chr(10),' | ')[-160:])
    return ok
ok = True
ok &= run('prove_S12.py')                 # S1,S4 lemmas, all (d,m)
ok &= run('prove_S34.py','3'); ok &= run('prove_S34.py','4')   # S''',S''''/S
ok &= run('prove_tops.py')                # top equations, all s
ok &= run('step_48.py')                   # V(3) non-adjacent
ok &= run('step_49.py')                   # V(4) non-adjacent (~6 min)
ok &= run('step_52.py')                   # slope (4,1) non-adjacent
for ds in [(3,0),(3,1),(4,0),(4,1),(4,2)]: ok &= run('prove_adj34.py',*ds)
ok &= run('prove_E3b34.py','3'); ok &= run('prove_E3b34.py','4')
print('ALL PASS' if ok else 'SOME FAIL')
# plus referee's own modules in this folder: tow2r_lib.py (first-principles Phi,
# Laurent extraction, T1); executed checks logged in progress.log:
#  S^l/S from genintegrand definition vs N_l/Pi_l: 172 cells d=4..12, 0 bad (l=1..4)
#  P_r Stirling tails: recurrence/anchors/values/divisibility, all True
#  Phi_p first principles vs binteg: key (4,0), delta 0..4, 5/5 exact
```

## 6. Files in `new_results/TOW2R_review/`

`progress.log` (per-step log), `tow2r_lib.py` (referee's first-principles Φ/Laurent/T1 machinery), `tow2_copy/` (TOW2 scripts+libs, read-only reruns), step_NN.py files produced by the review harness (reruns, transcript forensics, independent S/P_r/Φ checks). Evidence: every residual reported above was printed by my own runs; all were exact 0 (rational-function identities over Q(i)); the 172 S-derivative cells and 5 ρ cells were exact matches against definition resp. binteg data.