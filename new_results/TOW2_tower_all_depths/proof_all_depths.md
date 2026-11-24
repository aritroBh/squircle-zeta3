# TOW2 — proof document: the tower law at v = +i, depths δ ≤ 4 PROVED; mechanism at all depths identified

Status: every step below is either PROVED IN WRITING (all indices) or VERIFIED BY EXACT CODE
(symbolic rational identities → remainder 0, or exact QSI data checks; ranges stated).
This extends TOW1 (which proved δ ≤ 2) to δ ≤ 4, and develops the machinery toward all depths.

## 0. Setup (inherited from TOW1, all PROVED there)

Keys (p,q), p ≥ q ≥ 0, degree d = p+q ≥ 1, level k ≥ d.  Tower coefficient
u_{pq}^{(k)}(δ) = [(v−i)^{−(2d−2−δ)}] HigherPart(W_{pq}^{(k)}) at v = +i, depth δ ∈ {0,…,2d−4}.
Initial data ρ_{pq}(δ) = u_{pq}^{(d)}(δ).  Recursion T1 (PROVED in TOW1) expresses
u^{(k)}(δ) via degree-(d+1) parents at level k.  Newton formulation: the tower law at depth δ
(polynomiality in k of degree ⌊δ/2⌋) ⟺ the triangular family of identities
    (N(δ,s))   e(δ,s) + e(δ,s+1) = Λ^{(δ)}(s),     0 ≤ s ≤ ⌊δ/2⌋,
with e(δ,s) = 0 for s > ⌊δ/2⌋, e(δ,0) = ρ(δ), and Λ^{(δ)}(s) the T1 feed built from the
parents' e(δ′,s) (diagonal parents via the ê alternating conversion; fold for adjacent keys).
The diagonal key: ê_{pp}(δ,s) = (p+1)Σ_{δ′≤δ} τ1_{δ−δ′} e_{p+1,p}(δ′,s)/(4p−δ′−1)
(its parents are off-diagonal, so its law follows from the parent's; only the level k = 2p
is exceptional — handled by the blindness condition E3(δ)).

## 1. Initial data ρ(δ) at depths 3, 4 — PROVED (closed forms, all (d,m))

A(v) = Φ_p(v)(v−i)^{2d−2} analytic at i; ρ_{pq}(δ) = (1+[p>q])·A^{(δ)}(i)/δ!  (TOW1 Thm ρ).
logA = (2d−1−D)log v + D log(4v) + log S(t(v)) + (2−2d)log(v+i) − log(v²−a²) − log(v²−b²),
D = deg S = 2d, t(v) = (v²+2v−1)/(4v).

*Lemma T (t-derivatives, PROVED).*  t(i) = t* = (1+i)/2, t′(i) = 0,
t^{(r)}(i) = (−1)^{r−1} r!/(4 i^{r+1}) for r ≥ 2.  (t(v) = v/4 + 1/2 − 1/(4v): direct.)

*Lemma S1 (l = 1) and Lemma S4 (l = 2) — PROVED for all (d,m)* (`prove_S12.py`, remainders 0):
    S′(t*)/S(t*)  = N₁/(2d−3),          N₁ = 4m + (−2+6i)d − 4id²,
    S″(t*)/S(t*)  = N₂/((2d−3)(2d−5)),  N₂ = −40im + 16m² + (80/3+20i)d + (−16+96i)dm
                 + (−80−48i)d² − 32id²m + (208/3+16i)d³ − 16d⁴.
  Method: the normalized derivative T_l/T₀ is a finite sum (j from d−l to d) of explicit
  rational functions of (d,m) over Q(i): binom(1/2,j)/binom(1/2,d) and binom(j,s)/binom(d,m)
  are rational, w(s,m)/2^m = S(m,m−r)/m^underline(r) uses the Stirling-tail polynomials
  P_r(x) = S(x,x−r) (r ≤ 2: P₁ = x(x−1)/2, P₂ = x(x−1)(x−2)(3x−5)/24, PROVED for all x by
  the recurrence P_r(x+1) = (x+1−r)P_{r−1}(x) + P_r(x) verified symbolically + anchors), and
  (q2ⁿ)^{(a)}(t*), (B)^{(b)}(t*) are explicit (q2(t*) = 0, q2 = 2(t−t*)(t−t̄*)).
  sympy cancels the total minus N_l/Π_l to 0 — exact rational identity, all (d,m).
  (l = 1 reproves TOW1's Corollary S3: N₁/(2d−3) = −2id − 2(d−2m)/(2d−3).)

*Corollary (log-derivatives, PROVED).*  (log S)^{(l)}(t*) = sl_l:  sl₁ = R₁, sl₂ = R₂ − R₁²
with R_l = N_l/Π_l.  (log A)^{(r)}(i)/r! = L_r explicit: Faa di Bruno on log S∘t (only
blocks of size ≥ 2 since t′(i) = 0) + the algebraic pieces
(2d−1)(−1)^{r−1}i^{−r}/r + (2−2d)(−1)^{r−1}(2i)^{−r}/r − (−1)^{r−1}Σ_{c∈{±a,±b}}(i−c)^{−r}/r.

*Theorem ρ(3), ρ(4) (PROVED).*  ρ(δ)/(1+[p>q])·e00-factor = a_δ := [y^δ] exp(Σ L_r y^r)
(complete Bell polynomial):
    a₃ = (4id⁴ − 36id³ + (6+125i)d² − 12dm + (−18−198i)d + 36m + 117i)/(24d − 36),
    a₄ = (16d⁶ − 240d⁵ + (1528−48i)d⁴ + 96id³m + (−5296+432i)d³ − 864id²m + (10525−1308i)d²
          + (48+2616i)dm + (−11321+1320i)d − 48m² − 2640im + 5040)/(384d² − 1536d + 1440).
  (a₂ = TOW1's R₂, cross-checked symbolically.)  Exact check vs recorded initial towers:
  72 cells, all keys, d ≤ 11, residual 0 (`step_43`).

## 2. The zigzag top coefficients — PROVED (all s, given E0)

Data claim (TOW1 §7b's "Catalan" guess is wrong for s ≥ 4): the top Newton coefficients are
    e(2s,s)/e00   = (−1)^s A_{s+1} / Π_s,         Π_s := Π_{j=1}^s (2d−2j−1),
    e(2s+1,s)/e00 = (−1)^{s+1} i (d−s−2) A_{s+1} / Π_s  = −i(d−s−2)·e(2s,s)/e00,
where A_n = Euler zigzag numbers (sech x + tanh x = Σ A_n xⁿ/n!: 1,1,1,2,5,16,61,272,…).
  Exact check: s ≤ 6, all keys d ≤ 11 (`step_21`, `step_22`; note A₂=1, A₃=2, A₄=5 coincide
  with Catalan C₁,C₂,C₃ — that is why the wrong guess survived TOW1's data; A₅ = 16 ≠ 14).

*Theorem (top equations, PROVED for all s, non-adjacent keys; `prove_tops.py`).*
(i) Even top: e(2s,s) = Λ^{(2s)}(s).  Only the δ′ = 2s parent feeds contribute:
    Λ^{(2s)}(s) = B_d(2s)(e(2s,s)^parents) = [(p+1)τ1₀ rA₀ + (q+1)τ2₀ rB₀]·c_s(d+1)/(2d−2s−1)·e00
  with c_s(d) = (−1)^s A_{s+1}/Π_s.  The bracket equals (2d−1) — this is E0 (WT2b I1/I2,
  PROVED; `prove_tops.py` line 1).  The equation ⟺ Π_s(d+1)/Π_s(d) = (2d−1)/(2d−2s−1),
  exact telescoping (s−1 common factors cancel).  Note the numerator A_{s+1} cancels:
  the top equation fixes e(2s,s) only up to a (d,m)-independent constant; the zigzag
  constant itself is forced by the value equations V(2s) below.
(ii) Odd top: e(2s+1,s) = Λ^{(2s+1)}(s):  δ′ = 2s+1 and δ′ = 2s feeds.  With
  e(2s+1,s) = −i(d−s−2)c_s(d)e00 and M¹ := (p+1)τ1₁ rA₀ + (q+1)τ2₁ rB₀ = i(2d−1)(τ1₁−τ2₁)
  = (3i/2)(2d−1) (`prove_tops.py` line 2; τ1₁−τ2₁ = 3/2), the equation reduces (same
  telescoping) to  −(X−2)/(2X−1) = −(X−1)/(2X−2) + (3/2)/(2X−1), X = d−s,
  an exact identity (`prove_tops.py` line 3, remainder 0).
(iii) Adjacent keys: the diagonal-parent feed at Newton level s uses
  Σ_t(−1)^t ê(2s, s+t) = ê(2s,s) (only t = 0); the fold's level-s part at depth 2s
  likewise collapses; both top identities then hold: verified symbolically for all p at
  (δ,s) = (3,1), (4,2) (`prove_adj34.py`, remainders 0).  (At (2,1) this is TOW1's E2a.)

## 3. The Newton identities at depths 3, 4 — PROVED for all d

Closed forms of the full e-family at δ ≤ 4 (e(δ,0) = ρ(δ) of §1; tops of §2;
    e(4,1)/e00 = (d³ − 9d² + (113/4)d − id/2 + im − 29)/((2d−3)(2d−5))
found by exact fit over 32 (key,d) points; verified at fresh points d = 12,13).

*Theorem.*  For every key of degree d and both cases (non-adjacent p ≥ q+2; adjacent
p = q+1, with the diagonal ê-conversion and the fold), and every (δ,s) with δ ≤ 4:
    e(δ,s) + e(δ,s+1) = Λ^{(δ)}(s).
  Proof: exact rational identities in (d,p) resp. p over Q(i); remainder 0:
  - non-adjacent: V(3) = (3,0) [24 s], V(4) = (4,0) [281 s], E4-slope (4,1) [0.2 s]
    (`steps 51,52,55`; the tops by §2);
  - adjacent: all five (3,0),(3,1),(4,0),(4,1),(4,2) [0.3–51 s, `prove_adj34.py`].
  The τ-constants are all in Q(i) (τ1₀ = (−1−i)/2, τ1₁ = 3(1−i)/4, τ1₂ = 1+3i/4,
  τ1₃ = −5/8+5i/4, τ1₄ = −23/16−7i/16, τ2 = conj-√2-partner; verified, `step_58`).

*Theorem (blindness E3 at δ = 3, 4, PROVED for all p; `prove_E3b34.py`).*  The feed of
  Δ_pp(δ′) = ρ_pp(δ′) − Σ_s(−1)^s ê_pp(δ′,s) (δ′ ≤ δ) into the child (p,p−1) at depth δ
  (DL2-direct + fold) vanishes identically:
    Δ(3) = −(2p−3)(32p³−96p²+130p−69)/(12(4p−3)),
    Δ(4) = i(512p⁶−3840p⁵+12992p⁴−25600p³+30722p²−20789p+6054)/(48(4p−5)(4p−3)),
  and the feeds cancel to 0 [1.1 s, 2.8 s].

*Exact check of the whole tower law at δ ≤ 4:* closed-form candidate U vs the recorded
towers, all keys, k ≤ 14: 1348 cells, 0 bad (`step_77`).

## 4. Conclusion of the depth-4 proof

By TOW1's reduction theorem (PROVED there, §4): E0 (WT2b, PROVED) + the identities N(δ,s)
for δ ≤ 4 + blindness E3(δ) for δ ≤ 4 ⟹ **the tower law holds at depths δ ≤ 4 for all
keys and all k** (off-diagonal: k ≥ d; diagonal: k ≥ 2p+1 with the k = 2p exception):
u_{pq}^{(k)}(δ) is a polynomial in k of degree ⌊δ/2⌋ (exactly, by the zigzag tops).

## 5. Toward all depths (Task B status)

The mechanism is now fully explicit (see result.md for the master equation): the tower law
at all depths ⟺ the single generating identity  ρ(z) = Λ̂(Û(z,−1))  (value equations V(δ)
for all δ) + the top equations (PROVED for all s, §2) + blindness at all depths.
Components in place for δ = 5,6: R₃, R₄ = S‴/S, S⁗/S fitted and verified (all-(d,m)
symbolic proof pending in `prove_S34.py`); ρ(5), ρ(6) via the same Bell machinery
(L_r needs sl_l, l ≤ ⌊r/2⌋); the e(δ,s) closed forms are fitted as needed.  See result.md.


## Referee errata (TOW2R + CX8, both PROOF COMPLETE; recorded by the M1 session 2026-10-28; the text above is unchanged)
value/slope proofs are step_48/49/52 (not 51/52/55); R_3, R_4 are proved (prove_S34.py), not pending; normalise rho(delta) and Delta(delta) by rho(0) / rho_pp(0); separate off-diagonal e from stable diagonal e-hat; state the finite descending-induction uniqueness argument; zigzag = sec+tan; the all-depth top law and the master equation are outlook only (x^0 = initial condition, x^1 = value identity); tow2lib.e_desc does not terminate and step_47.py crashes (not load-bearing); domains: delta = 3, 4 need d >= 4, adjacent p >= 3; the tower->kernel map is real-linear.
