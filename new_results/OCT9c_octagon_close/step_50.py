# verify_oct9c.py -- OCT9c: numerical re-check of the P'''(1) closed form for the generic
# rational octagon K = {|x|+|y|+|2x+y|+|x+3y| <= 1}.  Run from OCT9c_octagon_close/ (< 3 min).
# (1) Evaluates P'''(1) = 12*sum_e |y'|_e*B_e at 240 digits FROM SCRATCH with mpmath
#     (polylog/log/pi/zeta only; no OCT9b value caches; the closed form has no Li21 atoms).
# (2) Compares with (a) the hardcoded 80-digit assembled target and
#     (b) independent per-edge quadrature of int_0^1 E3 dt at 110 digits.
import os, sys, pickle
import mpmath as mp, sympy as sp
from sympy import Rational as R
BASE = "<project>/research/08_theory"
O9B = os.path.join(BASE, "OCT9b_octagon_debug")
sys.path.insert(0, O9B); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e3build import quad_Earc
from oct9c_common import (load_all, get_engine, atom_values, rhs_terms_235, rhs_terms_2357,
                          eval_rhs_235, eval_rhs_2357)
DPS = 240; mp.mp.dps = DPS
z2 = mp.pi**2/6; z3 = mp.zeta(3); l2, l3, l5, l7 = mp.log(2), mp.log(3), mp.log(5), mp.log(7)
INV = {  # OCT9c task-1 invisible constants (PSLQ relations identical at 250 and 340 digits)
 0: (mp.mpf(30167467635454913)*z3 + mp.mpf(18544071891072555)*z2*l2 + mp.mpf(40489125837856272)*z2*l3
     + mp.mpf(15850884944620440)*z2*l5 - mp.mpf(2646857348685000)*z2*l7) / mp.mpf(654492976038000),
 1: (mp.mpf(113678924735323)*z3 + mp.mpf(19198429982352)*z2*l2 + mp.mpf(5597936751276)*z2*l3
     + mp.mpf(15502556024064)*z2*l5) / mp.mpf(3700404000000),
 2: (mp.mpf(1978761997948319)*z3 - mp.mpf(332709802532496)*z2*l2 + mp.mpf(211436534059812)*z2*l3
     + mp.mpf(1019740246328664)*z2*l5) / mp.mpf(41624512450560),
 3: (mp.mpf(15326623217074655)*z3 + mp.mpf(1322454819718512)*z2*l2 + mp.mpf(4740914439129636)*z2*l3
     + mp.mpf(9254478478368600)*z2*l5) / mp.mpf(142615945728000)}
LEN = {0: mp.sqrt(41)/20, 1: mp.mpf(5)/36, 2: 5*mp.sqrt(13)/72, 3: mp.sqrt(29)/40}
D = load_all(); Gargs, Dargs = rhs_terms_235(); sols7, negw7, logprods7 = rhs_terms_2357()
P3 = mp.mpf(0)
for e in range(4):
    eng = get_engine(e); NUM = atom_values(eng, e, D['lam'], DPS)  # lam_mine.pkl (360-digit Li21)
    if e == 0:
        Wlo = D['Wp0']['Wlo_e0']; RHS = eval_rhs_2357(D['s0']['solB0'], sols7, negw7, logprods7, DPS)
    else:
        Wlo = D['s123'][e]['Wlo']; RHS = eval_rhs_235(D['s123'][e]['sol'], Gargs, Dargs, DPS)
    B = mp.mpf(sp.N(Wlo.subs(list(NUM.items())), DPS)) + RHS + INV[e]
    P3 += 12*LEN[e]*B
    print(f"edge {e}: B_e = {mp.nstr(B, 30)}")
print("closed form P'''(1) =", mp.nstr(P3, 70))
TGT = mp.mpf('0.62257543954218437801774878116005493382629603171621955014299935560420360819318155')
print("|closed - assembled target| =", mp.nstr(abs(P3 - TGT), 3)); assert abs(P3 - TGT) < mp.mpf('1e-70')
EDGES = [
 (R(1,4),R(1,5),R(1,5),R(3,2),R(3,5),R(-7,12),[0,R(1,5)],[R(1,4),0]),
 (R(1,9),R(1,4),R(1,4),R(1,3),R(1,2),R(1,9),[R(1,4),0],[R(1,3),R(-1,9)]),
 (R(5,9),R(5,8),R(1,8),R(5,3),R(1,4),R(-5,9),[R(1,8),R(-1,4)],[R(1,3),R(-1,9)]),
 (R(1,8),R(1,5),R(1,5),R(1,4),R(3,5),R(1,24),[0,R(-1,5)],[R(1,8),R(-1,4)])]
mp.mp.dps = 110
Q = 12*sum(quad_Earc(*ed) for ed in EDGES)
print("quadrature P'''(1) (110 dps) =", mp.nstr(Q, 60)); mp.mp.dps = DPS
print("|closed - quadrature| =", mp.nstr(abs(P3 - Q), 3)); assert abs(P3 - Q) < mp.mpf('1e-90')
print("verify_oct9c: ALL CHECKS PASS")
