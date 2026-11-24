import sys, time, pickle
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import atom_values
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
LEN = {0: mp.sqrt(41)/20, 1: mp.mpf(5)/36, 2: 5*mp.sqrt(13)/72, 3: mp.sqrt(29)/40}
DPS = 300
mp.mp.dps = DPS
z2 = mp.pi**2/6; z3 = mp.zeta(3)
def inv_val(e):
    rel = INV[e]
    consts = [z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5)] + ([z2*mp.log(7)] if e==0 else [])
    return -mp.fsum([mp.mpf(rel[1+k])*consts[k] for k in range(len(consts))])/mp.mpf(rel[0])
P3 = mp.mpf(0)
Bvals = {}
for e in (0,1,2,3):
    eng = get_engine(e)
    NUM = atom_values(eng, e, lam_hi, DPS)
    if e == 0:
        Wlo = D['Wp0']['Wlo_e0']; c = sol0
        li = mp.fsum([mp.mpf(str(c[i]))*Li3v(B0[i]) for i in range(20)])
        pr = mp.fsum([mp.mpf(str(c[20+j]))*(-Li2v(PB0[j][0])*mp.log(PB0[j][1])) for j in range(44)])
    else:
        Wlo = D['s123'][e]['Wlo']; c = sols123[e]
        li = mp.fsum([mp.mpf(str(c[i]))*Li3v(G8[i]) for i in range(8)])
        pr = mp.fsum([mp.mpf(str(c[8+j]))*(-Li2v(PB19[j][0])*mp.log(PB19[j][1])) for j in range(19)])
    Wlo_v = mp.mpf(sp.N(Wlo.subs(list(NUM.items())), DPS))
    B = Wlo_v + li + pr + inv_val(e)
    Bvals[e] = (Wlo_v, li, pr, inv_val(e), B)
    P3 += 12*LEN[e]*B
    print(f"edge {e}: Wlo={mp.nstr(Wlo_v,12)} li={mp.nstr(li,12)} pr={mp.nstr(pr,12)} inv={mp.nstr(inv_val(e),12)} B={mp.nstr(B,18)} ({time.time()-t0:.1f}s)")
print("P'''(1) short form =", mp.nstr(P3, 70))
TGT = mp.mpf('0.62257543954218437801774878116005493382629603171621955014299935560420360819318155')
print("|short - OCT9c target| =", mp.nstr(abs(P3-TGT), 3))
pickle.dump({'P3_short_300': mp.nstr(P3, 300)}, open("P3_short.pkl","wb"))
