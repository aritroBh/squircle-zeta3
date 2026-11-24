import sys, time, pickle
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import atom_values, eval_rhs_235, eval_rhs_2357, rhs_terms_235, rhs_terms_2357
t0=time.time()
SF = pickle.load(open("short_form.pkl","rb"))
G8 = [sp.Rational(g) for g in SF['G8']]
PB19 = [(sp.Rational(w),p) for w,p in SF['PB19']]
B0 = [sp.Rational(g) for g in SF['B0']]
PB0 = [(sp.Rational(w),p) for w,p in SF['PB0']]
sols123 = {e: [sp.Rational(c) for c in SF['sols123'][e]] for e in (1,2,3)}
sol0 = [sp.Rational(c) for c in SF['sol0']]
INV = {e: pickle.load(open("invis_short.pkl","rb"))[e][1] for e in (0,1,2,3)}
lam_hi = pickle.load(open("lam_hi.pkl","rb"))
def Li3v(a): return mp.re(mp.polylog(3, mp.mpf(a)))
def Li2v(a): return mp.re(mp.polylog(2, mp.mpf(a)))
DPS = 300
mp.mp.dps = DPS
z2 = mp.pi**2/6; z3 = mp.zeta(3)
def inv_val(e):
    rel = INV[e]
    consts = [z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)] + ([z2*mp.log(7)] if e==0 else [])
    return -mp.fsum([mp.mpf(rel[1+k])*consts[k] for k in range(len(consts))])/mp.mpf(rel[0])
# OCT9c invisibles
l2,l3,l5,l7 = mp.log(2),mp.log(3),mp.log(5),mp.log(7)
INV9c = {
 0: (mp.mpf(30167467635454913)*z3 + mp.mpf(18544071891072555)*z2*l2 + mp.mpf(40489125837856272)*z2*l3
     + mp.mpf(15850884944620440)*z2*l5 - mp.mpf(2646857348685000)*z2*l7) / mp.mpf(654492976038000),
 1: (mp.mpf(113678924735323)*z3 + mp.mpf(19198429982352)*z2*l2 + mp.mpf(5597936751276)*z2*l3
     + mp.mpf(15502556024064)*z2*l5) / mp.mpf(3700404000000),
 2: (mp.mpf(1978761997948319)*z3 - mp.mpf(332709802532496)*z2*l2 + mp.mpf(211436534059812)*z2*l3
     + mp.mpf(1019740246328664)*z2*l5) / mp.mpf(41624512450560),
 3: (mp.mpf(15326623217074655)*z3 + mp.mpf(1322454819718512)*z2*l2 + mp.mpf(4740914439129636)*z2*l3
     + mp.mpf(9254478478368600)*z2*l5) / mp.mpf(142615945728000)}
Gargs235, Dargs235 = rhs_terms_235(); sols7, negw7, logprods7 = rhs_terms_2357()
for e in (0,1,2,3):
    eng = get_engine(e)
    NUM = atom_values(eng, e, lam_hi, DPS)
    if e == 0:
        Wlo = D['Wp0']['Wlo_e0']; c = sol0
        li = mp.fsum([mp.mpf(str(c[i]))*Li3v(B0[i]) for i in range(20)])
        pr = mp.fsum([mp.mpf(str(c[20+j]))*(-Li2v(PB0[j][0])*mp.log(PB0[j][1])) for j in range(44)])
        RHS9c = eval_rhs_2357(D['s0']['solB0'], sols7, negw7, logprods7, DPS)
    else:
        Wlo = D['s123'][e]['Wlo']; c = sols123[e]
        li = mp.fsum([mp.mpf(str(c[i]))*Li3v(G8[i]) for i in range(8)])
        pr = mp.fsum([mp.mpf(str(c[8+j]))*(-Li2v(PB19[j][0])*mp.log(PB19[j][1])) for j in range(19)])
        RHS9c = eval_rhs_235(D['s123'][e]['sol'], Gargs235, Dargs235, DPS)
    Wlo_v = mp.mpf(sp.N(Wlo.subs(list(NUM.items())), DPS))
    B_short = Wlo_v + li + pr + inv_val(e)
    B_9c = Wlo_v + RHS9c + INV9c[e]
    print(f"edge {e}: B_short = {mp.nstr(B_short, 30)}")
    print(f"         B_9c    = {mp.nstr(B_9c, 30)}")
    print(f"         |diff|  = {mp.nstr(abs(B_short-B_9c), 2)}  ({time.time()-t0:.1f}s)")
