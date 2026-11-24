# A3: MY OWN symbolic forced recursion in (d,m): base a_delta = first-principles series (a_mine.pkl), tau from MY DL1/DL2 series (core.TAU),
#     F(delta,s+1)=Lam(delta,s)-F(delta,s);  closing F(delta,floor(delta/2)+1)=0 for delta=0..8;  compare with tow4lib.forced and tow3lib.e.
import sys, time, pickle
sys.dont_write_bytecode=True
W='<project>/research/08_theory/TOW4_tower_depth8'
sys.path.insert(0, W)
import sympy as sp
import core
import tow2lib as T2, tow4lib as T4, tow3lib as T3
d, m = T2.d, T2.m; I = sp.I
def tosym(g): return sp.Rational(int(g.r.p),int(g.r.q)) + I*sp.Rational(int(g.i.p),int(g.i.q))
tau1 = [tosym(z) for z in core.TAU1]; tau2 = [tosym(z) for z in core.TAU2]
amine = {k: sp.sympify(v, locals={'d':d,'m':m,'I':I}) for k,v in pickle.load(open('a_mine.pkl','rb')).items()}
F = {}
def Fget(dp, s):
    if s > dp//2: return sp.Integer(0)
    return F[(dp,s)]
def Lam(delta, s):
    tot = sp.Integer(0)
    for dp in range(delta+1):
        fA = Fget(dp,s).subs({d:d+1, m:m+1}, simultaneous=True)
        fB = Fget(dp,s).subs({d:d+1, m:m}, simultaneous=True)
        tot += I*(2*d-1)/(2*d-dp-1) * (tau1[delta-dp]*fA - tau2[delta-dp]*fB)
    return tot
t0=time.time(); closing={}
for delta in range(0,9):
    F[(delta,0)] = sp.cancel(amine[delta])
    for s in range(0, delta//2+1):
        nxt = sp.cancel(sp.together(Lam(delta,s) - Fget(delta,s)))
        if s+1 <= delta//2:
            F[(delta,s+1)] = nxt
        else:
            closing[delta] = (nxt == 0)
            print(f'  my closing delta={delta}: F({delta},{s+1}) == 0 : {nxt==0}   ({time.time()-t0:.1f}s)', flush=True)
print('MY closing conditions delta<=8:', closing, flush=True)
pickle.dump({k:str(v) for k,v in F.items()}, open('F_mine.pkl','wb'))
# compare with tow4lib.forced (loads forced_F78.pkl) and tow3lib.e (delta<=6)
ok_forced = True
for (dl,s),val in sorted(F.items()):
    df = sp.cancel(sp.together(val - T4.forced(dl,s)))
    if df != 0: ok_forced=False; print('  DIFF vs tow4lib.forced', dl, s, df)
print('F_mine == tow4lib.forced (all delta<=8, all s):', ok_forced, f'({time.time()-t0:.1f}s)', flush=True)
ok3 = all(sp.cancel(sp.together(F[(dl,s)] - T3.e(dl,s)))==0 for dl in range(0,7) for s in range(0,dl//2+1))
print('F_mine == tow3lib.e (delta<=6):', ok3)
# tops == zigzag
okz = all(sp.cancel(sp.together(F[(dl,dl//2)] - T2.e_top(dl)))==0 for dl in range(0,9))
print('F_mine(top) == zigzag e_top (delta<=8):', okz)
