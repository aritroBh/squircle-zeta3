# verify_oct9d.py -- OCT9d: numerical re-check of the SHORT closed form of P'''(1) for the
# generic rational octagon K = {|x|+|y|+|2x+y|+|x+3y| <= 1}.  Run from OCT9d_octagon_short/ (< 5 min).
# (1) Evaluates P'''(1) = 12*sum_e |y'|_e * (Wlo_e + Lambda_e + Pi_e + eps_e) at 240 digits FROM
#     SCRATCH (mpmath polylog/log/pi/zeta only; no value caches; Wlo expressions loaded from OCT9b
#     storage are part of the formula).  Li3/product/invisible data hardcoded from short_form.pkl.
# (2) Compares with (a) the stored 340-digit quadrature value (quad330_Q1.pkl, 12*sum E_arc),
#     (b) a live recomputed 110-dps per-edge quadrature, (c) the OCT9c 80-digit target.
import os, sys, pickle
import mpmath as mp, sympy as sp
from sympy import Rational as R
BASE = "<project>/research/08_theory"
O9B = os.path.join(BASE, "OCT9b_octagon_debug")
sys.path.insert(0, O9B); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e3build import quad_Earc
from oct9d_common import load_all, get_engine
DPS = 240; mp.mp.dps = DPS
z2 = mp.pi**2/6; z3 = mp.zeta(3); l2, l3, l5, l7 = mp.log(2), mp.log(3), mp.log(5), mp.log(7)
def Li3v(a): return mp.re(mp.polylog(3, mp.mpf(a)))
def Li2v(a): return mp.re(mp.polylog(2, mp.mpf(a)))
G8  = [R(-9,16), R(4,9), R(1,5), R(1,3), R(-5,3), R(-1,5), R(5,4), R(2,5)]
PB19 = [(R(-1,3),2),(R(-1,3),3),(R(-1,3),5),(R(-1,4),2),(R(-1,4),3),(R(-1,4),5),(R(1,4),2),(R(1,4),3),
        (R(1,4),5),(R(-5,3),2),(R(-5,3),3),(R(-5,3),5),(R(-5,4),2),(R(-5,4),3),(R(-5,4),5),(R(1,5),3),
        (R(1,5),5),(R(5,8),2),(R(5,8),3)]
B0  = [R(-2,5),R(-3,2),R(-7,5),R(7,12),R(1,3),R(2,7),R(2,3),R(-3,5),R(5,6),R(1,5),
       R(1,7),R(6,7),R(3,7),R(3,5),R(7,15),R(-3,4),R(7,8),R(10,9),R(-7,8),R(-2,7)]
SF = pickle.load(open("short_form.pkl","rb"))
PB0 = [(R(w), p) for w, p in SF['PB0']]
sols123 = {e: [R(c) for c in SF['sols123'][e]] for e in (1,2,3)}
sol0 = [R(c) for c in SF['sol0']]
INV = {0: (37607500*z3 + 495501375*z2*l2 + 1398618*z2*l3 - 48773505*z2*l5 + 5972500*z2*l7, 25845375),
       1: (-100593252*z3 + 135960066*z2*l2 - 216979261*z2*l3 + 141579454*z2*l5, 45562500),
       2: (15564441600*z3 + 467608842*z2*l2 + 1991227897*z2*l3 - 8323927054*z2*l5, -16400517120),
       3: (2933760000*z3 + 8724711182*z2*l2 + 490424003*z2*l3 - 8331555450*z2*l5, 1248716800)}
LEN = {0: mp.sqrt(41)/20, 1: mp.mpf(5)/36, 2: 5*mp.sqrt(13)/72, 3: mp.sqrt(29)/40}
D = load_all()
P3 = mp.mpf(0)
for e in range(4):
    eng = get_engine(e)
    NUM = {}
    for s, key in eng.ATOMS.items():          # fresh mpmath atom values; Wlo has no Li21 atoms
        kind = key[0]
        if   kind == 'log':  NUM[s] = mp.log(mp.mpf(key[1]))
        elif kind == 'z2':   NUM[s] = z2
        elif kind == 'z3':   NUM[s] = z3
        elif kind == 'Li2':  NUM[s] = Li2v(key[1])
        elif kind == 'Li3':  NUM[s] = Li3v(key[1])
        elif kind == 'Li21': NUM[s] = mp.mpf(0)   # does not occur in Wlo
    if e == 0:
        Wlo = D['Wp0']['Wlo_e0']; c = sol0
        lam = mp.fsum([mp.mpf(str(c[i]))*Li3v(B0[i]) for i in range(20)])
        prd = mp.fsum([mp.mpf(str(c[20+j]))*(-Li2v(PB0[j][0])*mp.log(PB0[j][1])) for j in range(44)])
    else:
        Wlo = D['s123'][e]['Wlo']; c = sols123[e]
        lam = mp.fsum([mp.mpf(str(c[i]))*Li3v(G8[i]) for i in range(8)])
        prd = mp.fsum([mp.mpf(str(c[8+j]))*(-Li2v(PB19[j][0])*mp.log(PB19[j][1])) for j in range(19)])
    num, den = INV[e]
    B = mp.mpf(sp.N(Wlo.subs(list(NUM.items())), DPS)) + lam + prd - num/mp.mpf(den)
    P3 += 12*LEN[e]*B
    print(f"edge {e}: B_e = {mp.nstr(B, 30)}")
print("short closed form P'''(1) =", mp.nstr(P3, 70))
# (a) stored 340-dps quadrature (12*sum E_arc; length factors NOT applied -- B_e already contains 1/|y'|)
Q340 = 12*sum(mp.mpf(s) for s in pickle.load(open("quad330_Q1.pkl","rb")))
print("|closed - stored 340dps quadrature| =", mp.nstr(abs(P3 - Q340), 3)); assert abs(P3 - Q340) < mp.mpf('1e-230')
# (b) live 110-dps quadrature
EDGES = [
 (R(1,4),R(1,5),R(1,5),R(3,2),R(3,5),R(-7,12),[0,R(1,5)],[R(1,4),0]),
 (R(1,9),R(1,4),R(1,4),R(1,3),R(1,2),R(1,9),[R(1,4),0],[R(1,3),R(-1,9)]),
 (R(5,9),R(5,8),R(1,8),R(5,3),R(1,4),R(-5,9),[R(1,8),R(-1,4)],[R(1,3),R(-1,9)]),
 (R(1,8),R(1,5),R(1,5),R(1,4),R(3,5),R(1,24),[0,R(-1,5)],[R(1,8),R(-1,4)])]
mp.mp.dps = 110
Q = 12*sum(quad_Earc(*ed) for ed in EDGES)
mp.mp.dps = DPS
print("|closed - live 110dps quadrature| =", mp.nstr(abs(P3 - Q), 3)); assert abs(P3 - Q) < mp.mpf('1e-90')
# (c) OCT9c target
TGT = mp.mpf('0.62257543954218437801774878116005493382629603171621955014299935560420360819318155')
print("|closed - OCT9c target| =", mp.nstr(abs(P3 - TGT), 3)); assert abs(P3 - TGT) < mp.mpf('1e-70')
print("verify_oct9d: ALL CHECKS PASS")
