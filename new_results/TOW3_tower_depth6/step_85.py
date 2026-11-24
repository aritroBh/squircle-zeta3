work = '<project>/research/08_theory/TOW3_tower_depth6'
with open(f'{work}/progress.log','a') as f:
    f.write('step_81: Task D: -a map VALIDATED vs SUB1 table d=2..5 (0 fails); M_d from map d=2..7 ALL OK; S_d table computed; FINAL written\n')
# append Task D section to proof_depth6.md
proof_add = r"""
## 5. Task D — the sub-top map and M_d (DERIVED, validated; M_d closed form NOT found)

*The map (derived, validated).* For every key (i,j) of degree d (write n = d+1),
    Res_{v=-a} K^{(k)}_{ij} = -(i+1) Gm_{i+1,j}(-a)
(Res_{-a} DL1 = 1, Res_{-a} DL2 = 0, the fold -iF(Gm_ii DL1) has no pole at -a: its poles are
at 1/{poles of Gm_ii DL1} subset {+-i, a, -b}). With e_m = -h(m+1)/m the integrated tower
coefficients and w_m = (-a-i)^{-m} - (b-i)^{-m} (always sqrt2 x Gaussian rational),
    Gm(-a) = 2 Re_i[ sum_m e_m w_m ],  h^{(k)}(r) = e00 sum_s e(2n-2-r, s) binom(k-n, s).
Since beta_{d,M}(r) = Delta^{M-1} c_r(d+1) picks out the binom(k-n, M-1) coefficient:
    beta_{d,M}(-a) = -(i+1) * 2 Re_i[ e00 * sum_{m: floor((2d-m-1)/2) >= M-1} (-1/m) e(2d-m-1, M-1) w_m ].
Only the sqrt2-part Es of e00 contributes to the rational residue, and it cancels in the ratio.
VALIDATED against SUB1's exact table (top AND sub-top at -a, all keys, d = 2..5, 0 failures;
kernel-level check: the beta's reproduce the recorded kern_k residues, 4/4 spot checks).

*M_d from the map.* With S_d := e(2d-2, d-2) at (d+1, d+1) = (s0, s1) (Gaussian rationals),
o1 = (-1)^{d-1} A_{d-1}/Pi_{d-2}(d+1), e0 = (-1)^{d-2} A_{d-1}/Pi_{d-2}(d+1), n = d+1:
    M_d = 4 (-1)^{d-1} (2d-1)!! * B_d,
    B_d = (s0 - s1)/2 - o1/4 - e0/6   (d even);   B_d = (s0 + s1)/2 - o1/4 - e0/12  (d odd).
VALIDATED: d = 2..7 gives M_d = -1, 10, 15, 131, 394, 3151 — exactly SUB1's values (the
d = 6, 7 inputs S_6 = e(10,4)@(7,7), S_7 = e(12,5)@(8,8) are exact from the WTS towers).
S_d table (exact): d=2: (-1/2,-1/2); d=3: (4/15,2/15); d=4: (-9/70,-1/21); d=5: (37/630,1/63);
d=6: (-29/1155,-8/1485); d=7: (557/54054,244/135135). This is a genuine REDUCTION of M_d to a
single sub-top tower value per d (previously M_d was an unanalysed intercept).

*The sub-top family.* The closed forms (s = 2,3,4: e(4,1) PROVED TOW2, e(6,2) PROVED here,
e(8,3) CANDIDATE verified 411 fresh keys 0 bad) share the rigid shape
    e(2s, s-1) = [ (-1)^s C_{s-1} (d^3 + i m) + a_s d^2 + b_s d + c_s ] / Pi_s,
C_{s-1} Catalan (1, 2, 5): leading d^3 and i m coefficients coincide and are Catalan — a new
structural observation. e(10,4), e(12,5) are extractable from WTS towers (not fitted here).

*M_d closed form.* NOT found (consistent with SUB1's extensive negative search: no 3-term
zigzag recurrence — tested here too: alpha A_{d+1} + beta A_d + gamma A_{d-1} fails at d = 5).
The map above reduces the problem to the sub-top family a_s, b_s, c_s (rational-in-s data
s = 2..4 known, s = 5,6 extractable); closing a_s, b_s, c_s in s is the remaining step.
Status: CONJECTURE-level for the family shape; the map formula itself is exact data
(validated d <= 7 against SUB1's d <= 10 table where it overlaps).
"""
with open(f'{work}/proof_depth6.md','a') as f: f.write(proof_add)
print('done')
