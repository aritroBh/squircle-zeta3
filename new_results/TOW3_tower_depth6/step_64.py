work = '<project>/research/08_theory/TOW3_tower_depth6'
with open(f'{work}/progress.log','a') as f:
    f.write('step_61: Task C: degree-4 kernels k=5 vs cubic(k=6..9) 40/40; cubicity k=6..13 40/40\n')
result = r"""# TOW3 — the tower law at v = i: **PROVED TO DEPTH 6** (all d, all k); weight-5 polynomial law PROVED for every k >= 5; sub-top/M_d map (Task D)

## VERDICT

**PROVED TO DEPTH 6.** The tower law "the depth-delta tower coefficient u_{pq}^{(k)}(delta) at v = +i is a
polynomial in k of degree floor(delta/2)" is proved for depths **delta <= 6**, all keys (p,q) of all degrees d,
all k (off-diagonal: k >= d; diagonal: k >= 2p+1, k = 2p exceptional), extending TOW2 (delta <= 4) by exactly
this job's request: the value identities **V(5) and V(6)** — TOW2's stated sole gap — are here PROVED as exact
rational identities in (d,m) (numerator 0), together with the slopes, tops, adjacent-key identities (fold +
ehat-conversion) and the blindness E3 at depths 5 and 6.

**Consequence (Task C): the weight-five polynomial law for every k.** Tower law at degree 5 to depth
6 = 2*5-4 (proved here) ==(WTS Lemma 3, PROVED)==> kernel law K(4): degree-4 kernel residues are polynomials
in k of degree <= 3 for k >= 6 ==(WTS Lemma 2, k-independent linear image)==> W_5(k) is a cubic in k for
k >= 6. The level k = 5 is outside Lemma 3's range; it is here checked EXACTLY: the recorded degree-4 kernels
at k = 5 equal the cubic fitted from k = 6..9 at every pole, every QSI component, all three keys
(**40/40 components, 0 bad**; cubicity on k = 6..13 independently confirmed, 40/40). Hence
**W_5(k) = W5'' + (k-5) E1 + C(k-5,2) E2 + C(k-5,3) E3 for every k >= 5** (k = 4 exceptional, as always;
constants E1, E2, E3, W5'' as determined by WT5). This upgrades WT5's exact-data law (5 <= k <= 14) to a
proved law for all k >= 5.

**Task D (sub-top coefficients and M_d):** see Section D below (in progress).

## What is proved — every step in writing (all indices) or exact symbolic code (remainder 0)

| item | statement | status |
|---|---|---|
| rho(5), rho(6) closed forms | a_5 = N_5/(480(2d-5)(2d-3)), a_6 = N_6/(5760(2d-7)(2d-5)(2d-3)), Bell assembly from PROVED R_1..R_4 | PROVED + exact check **720 cells each, 0 bad** |
| V(5), V(6) non-adjacent | rho(delta)+e(delta,1) = Lambda^{(delta)}(0), exact identity in (d,m) | PROVED, numerator 0 (0.8s, 7.3s) |
| slopes (5,1), (6,1), (6,2) non-adjacent | e(delta,s)+e(delta,s+1) = Lambda^{(delta)}(s) | PROVED, numerator 0 (0.2-0.6s) |
| tops (5,2), (6,3) non-adjacent | e(delta,top) = Lambda^{(delta)}(top); zigzag | PROVED (here 0.0s; also TOW2 prove_tops all s) |
| adjacent identities (5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3) | with ehat-conversion + fold, rational identities in p | PROVED, remainder 0 (0.8-668 s, prove_adj56.py); (5,2),(6,3) are NEW adjacent tops s=2,3 |
| blindness E3 at delta = 5, 6 | Delta_pp feed into child (p,p-1) vanishes identically | PROVED all p, numerator 0 (9.4s, 20.9s, prove_E3b56.py) |
| e(5,1), e(6,1), e(6,2) closed forms | fitted (e(6,1) here; e(5,1),e(6,2) TOW2 fits) then PROVED via every identity they enter + 414 fresh data keys each | PROVED |
| tower law delta <= 6 | above + TOW1 reduction theorem | PROVED; exact check **38939 cells, 0 bad** (all keys d <= 52, k <= d+8) |
| weight-5 law (Task C) | K(4) for k >= 6 (proved) + k = 5 exact check 40/40 | PROVED for all k >= 5 |

## Key structures (all inherited from TOW2, two new closed forms)

- a_5, a_6: proof_depth6.md Section 1 (explicit numerators N_5, N_6 over Q(i)).
- e(6,1) = N_{61}/(96 Pi_3): proof_depth6.md Section 2; e(5,1) = N_{51}/(12 Pi_2), e(6,2) = N_{62}/(2 Pi_3).
- The value identities V(5), V(6) close: the sole gap of TOW2 is closed at depths 5, 6.

## Files (research/08_theory/TOW3_tower_depth6/)

- `tow3lib.py` (e-family delta <= 6, Lambda feed, arho via Bell from L_parsed.pkl),
  `prove_V56.py` (non-adjacent identities), `prove_adj56.py` (adjacent), `prove_E3b56.py` (blindness),
  `L_parsed.pkl`, `arho56_simplified.pkl`, `e61_closed.pkl`.
- `proof_depth6.md`, `result.md` (this file), `progress.log`. `verify_tow3.py` (in progress).

## Honest boundaries (not proved)

- Value identities V(delta) for delta >= 7 (proved <= 6); tower law at depths >= 7.
- Adjacent-key top equations for s >= 4 (proved now for (3,1),(4,2),(5,2),(6,3); open s >= 4).
- Task D's M_d closed form: in progress (candidates to be checked against SUB1 d = 2..10).
"""
with open(f'{work}/result.md','w') as f: f.write(result)
print('result.md updated')
