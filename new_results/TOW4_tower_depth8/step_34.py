work = '<project>/research/08_theory/TOW4_tower_depth8'
result = """# TOW4 — the tower law at v = i at depths 7 and 8; the weight-six polynomial law; all-depth attempt

## VERDICT (updated as steps finish)

**TASK A COMPLETE (rho(7), rho(8) PROVED).** Depths 7, 8 Newton identities + E3: IN PROGRESS.

## Task A — rho(7), rho(8) closed forms, PROVED

Ingredients: a_7 = Bell_7(L_1..L_7), a_8 = Bell_8(L_1..L_8) use sl_l, l <= 4 (R_1..R_4 PROVED in TOW2
`prove_S12.py`, `prove_S34.py`; rerun green by TOW3R). The simplified cache L(1..8) (`L_simplified.pkl`,
inherited from TOW2 via TOW3) is here verified EQUAL to the proved L(r) construction of tow2lib
**symbolically** (exact difference == 0) for **r = 1..8** — closing the r = 7, 8 gap left at TOW3
(TOW3R's step_04 did r <= 6). Hence rho(7), rho(8) rest entirely on proved ingredients.

  a_7 = N_7 / (40320 (2d-7)(2d-5)(2d-3))     [N_7: Gaussian-integer polynomial, deg 10 in d, 41 terms]
  a_8 = N_8 / (645120 (2d-9)(2d-7)(2d-5)(2d-3))   [N_8: deg 12 in d, 57 terms]
  (Pi-part of denominators: Pi_3, Pi_4 = prod_{j=1}^{l}(2d-2j-1), matching a_delta ~ Pi_{floor(delta/2)}.)

Exact check vs TOWN1 towers_D8_K52.pkl: rho(delta)/rho(0) = a_delta with rho(0) = e00 (the 1+[p>q]
factor is inside e00; the initial level is k0 = d for ALL keys — for diagonal keys this is the
exceptional level k = 2p, whose rho_pp is the true initial data; sqrt2-part of ratio asserted 0):
**717 cells each for delta = 7, 8, all keys with 2d-4 >= delta, d <= 52; 0 bad.** (Sanity reproduces
TOW3's a_5, a_6 720-cell checks.)

## Task B — Newton identities at depths 7, 8 (IN PROGRESS)

Strategy (per brief): the FORCED ascending recursion F(delta,0) := a_delta, F(delta,s+1) :=
Lambda(delta,s) - F(delta,s) (fit-free, from proved rho and tau constants); the tower law at depth
delta is equivalent to the closing condition F(delta, floor(delta/2)+1) = 0 together with E3
(TOW3 claude_screen M1 pre-screen). To prove: F(7,4) = 0, F(8,5) = 0 symbolically in (d,m);
adjacent-key identities (p symbolic, fold + ehat-conversion); E3 at delta = 7, 8.

## Task C — weight-six law (pending Task B)

## Task D — all-depth attempt (pending)
"""
open(work + '/result.md','w').write(result)
proof = """# TOW4 — proof document: the tower law at depths 7 and 8

Inherited: everything TOW2 (delta <= 4) and TOW3 (delta <= 6) proved, including R_1..R_4,
the Bell theorem for rho, the zigzag tops, the TOW1 reduction (E0 + N(delta,s) + E3 ==> tower law),
the tau constants in Q(i), the adjacent machinery (fold + ehat; e00_pp anti-real so Dc = -conj(Delta);
e00_{p,p-1} real), and the FORCED-recursion equivalence (TOW3 claude_screen): the tower law at depth
delta <=> closing condition F(delta, floor(delta/2)+1) = 0 + E3, with F the fit-free ascending family
from rho.

Domains (carried forward, per TOW3R wording fix): depth delta needs 2d-4 >= delta (delta = 7, 8
=> d >= 6); adjacent (p,p-1) needs 4p-6 >= delta (p >= 3); k >= d off-diagonal; k >= 2p+1 diagonal
stable, k = 2p exceptional.

## 1. rho(7), rho(8) — PROVED (Task A)

a_7 = Bell_7, a_8 = Bell_8 from L(1..8); the L cache equals tow2lib.L(r) symbolically, r = 1..8
(diff == 0; r = 8: 9.9 s). Exact check: 717 cells each, 0 bad (k0 = d all keys; a_delta in Q(i)).

## 2. The Newton identities at depths 7, 8 — via the forced recursion (IN PROGRESS)
"""
open(work + '/proof_depth8.md','w').write(proof)
print("wrote result.md and proof_depth8.md")
