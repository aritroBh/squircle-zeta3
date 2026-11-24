work = '<project>/research/08_theory/TOW3_tower_depth6'
result = r"""# TOW3 — the tower law at v = i: depths 5 and 6, the weight-five law, and the sub-top/M_d map

## VERDICT (updated as the job progresses; currently PARTIAL)

**Task A COMPLETE.** rho(5), rho(6) assembled in closed form from PROVED ingredients
(the S-derivative lemmas R_1..R_4 of TOW2 `prove_S12.py`/`prove_S34.py` and the proved
Bell/Faa-di-Bruno theorem rho(delta)/e00 = a_delta = [y^delta] exp(sum L_r y^r)) and
VERIFIED against the exact TOWN1 towers: **720 cells each, 0 bad** (all keys, degrees
d <= 52 with 2d-4 >= delta; sanity: TOW2-proved rho(4) also 723/723).

  a_5 = N_5/(480 (2d-5)(2d-3)),  a_6 = N_6/(5760 (2d-7)(2d-5)(2d-3))
with N_5, N_6 explicit polynomials in (d, m) over Q(i) (proof_depth6.md, Section 1).

**Tasks B–E in progress.**

## Status table
| item | status |
|---|---|
| a_5, a_6 closed forms (Bell assembly from proved R_1..R_4) | PROVED + 720 exact cells each |
| value identities V(5), V(6) | in progress |
| slopes (5,1),(6,1),(6,2), tops (5,2),(6,3) | tops already PROVED all s (TOW2 prove_tops); slopes in progress |
| adjacent keys at delta 5,6 | in progress |
| blindness E3 at delta 5,6 | in progress |
| weight-five law (Task C) | pending B |
| sub-top e(2s,s-1), M_d (Task D) | pending |
"""
with open(f'{work}/result.md','w') as f: f.write(result)
proof = r"""# TOW3 — proof document: depths 5 and 6

Inherited from TOW2 (all PROVED there): setup, Lemma T, Lemmas S1/S4/S3/S5
(R_1..R_4 = S^{(l)}(t*)/S(t*) = N_l/Pi_l for all (d,m)), the Bell theorem
rho(delta)/e00 = a_delta = [y^delta] exp(sum_{r>=1} L_r y^r) with L_r explicit from
sl_l (l <= floor(r/2)) and algebraic pieces, zigzag tops, the master reduction:
tower law at depth delta  <=>  { N(delta,s): e(delta,s)+e(delta,s+1) = Lambda^{(delta)}(s) }
+ blindness E3(delta); tops N(delta, floor(delta/2)) PROVED all s (non-adjacent).

## 1. rho(5), rho(6) — PROVED closed forms

Ingredients: a_5 = Bell_5(L_1..L_5) uses sl_l, l <= 2 (R_1, R_2 PROVED);
a_6 = Bell_6(L_1..L_6) uses sl_l, l <= 3 (R_1, R_2, R_3 PROVED).
L(1..8) simplified cache verified equal to the proved L(r) construction
(random-point exact check, r = 5, 6; the L(r) formula itself is the proved
Faa di Bruno assembly of TOW2 proof_all_depths.md Section 1).

  a_5 = ( -16 i d^7 + 304 i d^6 + (-80-2520 i) d^5 + 160 d^4 m + (1040+11840 i) d^4
          - 2080 d^3 m + (-5140-34069 i) d^3 + (10280-240 i) d^2 m + (11440+60081 i) d^2
          + 240 i d m^2 + (-22880+960 i) d m + (-9600-59920 i) d - 960 i m^2 + 19200 m + 25200 i )
        / (480 (2d-5)(2d-3))
  a_6 = ( -64 d^9 + 1728 d^8 + (-20976+480 i) d^7 - 960 i d^6 m + (150432-10080 i) d^6
          + 20160 i d^5 m + (-702876+89040 i) d^5 + (-2880-178080 i) d^4 m + (2221452-423840 i) d^4
          + 2880 d^3 m^2 + (34560+847680 i) d^3 m + (-4760069+1146870 i) d^3
          - 34560 d^2 m^2 + (-138960-2293260 i) d^2 m + (6677868-1670730 i) d^2
          + (138960-1440 i) d m^2 + (186480+3341460 i) d m + (-5524995+1019460 i) d
          + 960 i m^3 - 186480 m^2 - 2038920 i m + 1984500 )
        / (5760 (2d-7)(2d-5)(2d-3))

Exact check vs TOWN1 towers_D8_K52.pkl: rho(delta)/e00 = a_delta, all keys with
2d-4 >= delta, d <= 52: rho(5) 720/720, rho(6) 720/720, rho(4) (sanity) 723/723.
(Method: QSI division of the recorded values; sqrt2-part of the ratio = 0 asserted.)

## 2. The Newton identities at depths 5, 6 (in progress)
"""
with open(f'{work}/proof_depth6.md','w') as f: f.write(proof)
with open(f'{work}/progress.log','a') as f:
    f.write('step_33: rho(5),rho(6) VERIFIED 720 cells each 0 bad (Task A done); wrote result.md/proof_depth6.md v1\n')
print('written')
