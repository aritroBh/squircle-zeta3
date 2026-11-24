"""verify_tow4r.py -- TOW4R referee verifier (<10 min). Run with project venv python in TOW4R_review/.
PASS/FAIL per item. Exact arithmetic throughout; residual 0 means exactly 0.
Fast subset of this review's checks; slow certificates (adj 7,0/8,0 ~50-90s each, from-scratch
L_r 88 keys ~20 min) were run and logged separately -- see result.md and progress.log.
"""
import os, sys, subprocess, pickle, time
from fractions import Fraction as Fr
from math import comb
W = os.path.dirname(os.path.abspath(__file__))
TC = W + "/tow4_copy"
sys.path.insert(0, TC)
PY = sys.executable
RES = []
def rep(name, ok, extra=""):
    RES.append((name, bool(ok))); print(f'[{"PASS" if ok else "FAIL"}] {name}  {extra}', flush=True)
import sympy as sp
import tow2lib as T2, tow4lib as T4
d_, m_ = T2.d, T2.m
i = sp.I
t00 = time.time()
# 1. TAU8 from the DL definition (own exact computation)
aq = sp.sqrt(2)+1; bq = sp.sqrt(2)-1
t1 = (-1)**9*((i-bq)**(-10)+(i+aq)**(-10)-i**(-10)); t2 = (-1)**9*((i-aq)**(-10)+(i+bq)**(-10)-i**(-10))
rep("TAU8", sp.simplify(t1-T4.tau1f(8))==0 and sp.simplify(t2-T4.tau2f(8))==0, "exact vs DL def")
# 2. L cache == tow2lib.L(r), r=1..8, symbolic
Ls = {k: sp.sympify(v, locals={'d':d_,'m':m_,'I':i}) for k,v in pickle.load(open(TC+"/L_simplified.pkl",'rb')).items()}
ok = all(sp.cancel(sp.together(T2.L(r)-Ls[r]))==0 for r in range(1,9))
rep("L r=1..8 symbolic", ok)
# 3. forced recursion + closing F(7,4)=0, F(8,5)=0 (rerun of TOW4 certificate, in tow4_copy)
r = subprocess.run([PY,"-B","forced78.py"], capture_output=True, text=True, cwd=TC, timeout=600)
rep("CLOSE F(7,4)=0,F(8,5)=0", "VERDICT non-adjacent depth<=8: True" in r.stdout)
# 4. adjacent fast certificates: the two tops + (8,3) (all 9 rerun green in review; these cover both parities+tops)
ok = True
for dl,s in [(7,3),(8,4),(8,3)]:
    rr = subprocess.run([PY,"-B","prove_adj78_fast.py",str(dl),str(s)], capture_output=True, text=True, cwd=TC, timeout=600)
    ok &= (rr.returncode==0 and "== 0 : True" in rr.stdout)
rep("ADJ (7,3),(8,4),(8,3)", ok, "all 9 rerun green separately")
# 5. E3 certificates dl=7,8
ok = True
for dl in [7,8]:
    rr = subprocess.run([PY,"-B","prove_E3b78.py",str(dl)], capture_output=True, text=True, cwd=TC, timeout=600)
    ok &= (rr.returncode==0 and "== 0 : True" in rr.stdout)
rep("E3 dl=7,8", ok)
# 6. fast-method non-vacuity: perturbed feed must FAIL
src = open(TC+"/prove_adj78_fast.py").read()
p1 = src.replace("terms.append(p*tau2(delta-dp)*Ediag(dp, s)/(2*d-dp-1))",
                 "terms.append(-p*tau2(delta-dp)*Ediag(dp, s)/(2*d-dp-1))")
open(TC+"/_perturb_check.py","w").write(p1)
rr = subprocess.run([PY,"-B","_perturb_check.py","7","2"], capture_output=True, text=True, cwd=TC, timeout=600)
rep("non-vacuity (sign-flipped feed fails)", rr.returncode!=0 and "== 0 : False" in rr.stdout)
# 7. full-law data check (49953+222 cells)
rr = subprocess.run([PY,"-B","check_full78.py"], capture_output=True, text=True, cwd=TC, timeout=900)
out = rr.stdout
rep("LAW 49953+222 cells", "cells = 49953" in out and out.count("bad = 0")>=2)
# 8. RHO 717 cells each dl=7,8 with referee's own Q(r2,i) division
H = pickle.load(open("<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl",'rb')); H.pop('_meta',None)
def mydiv(x,y):
    def mul(A,B): return (A[0]*B[0]+2*A[1]*B[1], A[0]*B[1]+A[1]*B[0])
    Y=((y[0],y[1]),(y[2],y[3])); Yc=((Y[0][0],-Y[0][1]),(Y[1][0],-Y[1][1]))
    Nre=(mul(Y[0],Yc[0])[0]-mul(Y[1],Yc[1])[0],0); Nim=(mul(Y[0],Yc[1])[0]+mul(Y[1],Yc[0])[0],0)
    den=Nre[0]**2+Nim[0]**2; assert den!=0
    X=((x[0],x[1]),(x[2],x[3]))
    def cmul(A,B):
        return ((mul(A[0],B[0])[0]-mul(A[1],B[1])[0], mul(A[0],B[0])[1]-mul(A[1],B[1])[1]),
                (mul(A[0],B[1])[0]+mul(A[1],B[0])[0], mul(A[0],B[1])[1]+mul(A[1],B[0])[1]))
    Wd = cmul(cmul(X,Yc), ((Nre[0]/den,Fr(0)),(-Nim[0]/den,Fr(0))))
    return Wd
bad=0; n=0
for dl in [7,8]:
    aa = T4.arho(dl)
    for (p,q) in H:
        dd = p+q
        if 2*dd-4 < dl: continue
        Wd = mydiv(H[(p,q)][dl][dd], H[(p,q)][0][dd]); n+=1
        cv = sp.expand(aa.subs({d_:dd,m_:p}))
        if Wd[0][1]!=0 or Wd[1][1]!=0 or Wd[0][0]!=Fr(sp.re(cv)) or Wd[1][0]!=Fr(sp.im(cv)): bad+=1
rep("RHO 2x717 cells own division", bad==0, f"{n} cells, bad={bad}")
# 9. W6 kernel checks (own code, WTS data)
wts = "<project>/research/08_theory/WTS_weight_k_structure/theory"
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(6,14)}
keys5 = [(5,0),(4,1),(3,2)]
def lag(ys,xs,x):
    out = Fr(0)
    for a,xa in enumerate(xs):
        La = Fr(1)
        for b,xb in enumerate(xs):
            if b!=a: La *= Fr(x-xb, xa-xb)
        out += ys[a]*La
    return out
w6a = all(lag([K[k][key]['real'].get(pole,(Fr(0),)*4)[c] for k in [7,8,9,10,11]],[7,8,9,10,11],6)
          == K[6][key]['real'].get(pole,(Fr(0),)*4)[c]
          for key in keys5 for pole in set().union(*[set(K[k][key]['real']) for k in range(6,12)]) for c in range(4))
w6b = True
for key in keys5:
    for pole in set().union(*[set(K[k][key]['real']) for k in range(7,14)]):
        for c in range(4):
            cur = [K[k][key]['real'].get(pole,(Fr(0),)*4)[c] for k in range(7,14)]
            for _ in range(5): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
            if any(cur): w6b=False
w6c = all(not K[k][key]['ipos'] and not K[k][key]['ineg'] for k in range(6,14) for key in keys5)
rep("W6 k=6 on quartic(7..11) 48/48; quarticity 7..13; +-i poles empty 6..13", w6a and w6b and w6c)
# 10. MASTER telescoping at dl=7 (feed(U)-a_7==0 and identity with perturbed seed at dl=3,5)
eps = sp.Symbol('eps')
tau1 = T4.tau1f; tau2 = T4.tau2f
def build(delta, pert):
    Fz = {}
    def F(dp,s):
        return sp.Integer(0) if s>dp//2 else Fz[(dp,s)]
    def Lam(dl,s):
        tot = sp.Integer(0)
        for dp in range(0,dl+1):
            fA = F(dp,s).subs({d_:d_+1,m_:m_+1},simultaneous=True); fB = F(dp,s).subs({d_:d_+1,m_:m_},simultaneous=True)
            tot += i*(2*d_-1)/(2*d_-dp-1)*(tau1(dl-dp)*fA - tau2(dl-dp)*fB)
        return tot
    for dl in range(0,delta+1):
        Fz[(dl,0)] = sp.cancel(T4.arho(dl) + (eps if (pert and dl==delta) else 0))
        for s in range(0,dl//2): Fz[(dl,s+1)] = sp.cancel(Lam(dl,s)-F(dl,s))
        Fz[(dl,dl//2+1)] = sp.together(Lam(dl,dl//2)-F(dl,dl//2))
    return Fz
Fz = build(7, False)
okm = sp.expand(Fz[(7,4)])==0
Fz3 = build(3, True)
top=1
U = {dp: sum((-1)**s*(Fz3[(dp,s)] if s<=dp//2 else sp.Integer(0)) for s in range(0,dp//2+1)) for dp in range(0,4)}
feedU = sum(i*(2*d_-1)/(2*d_-dp-1)*(tau1(3-dp)*U[dp].subs({d_:d_+1,m_:m_+1},simultaneous=True)
        - tau2(3-dp)*U[dp].subs({d_:d_+1,m_:m_},simultaneous=True)) for dp in range(0,4))
tel = sp.expand(sp.together(Fz3[(3,2)] - (-1)**top*(feedU-(T4.arho(3)+eps))))
sens = sp.cancel(sp.diff(Fz3[(3,2)], eps)) != 0
rep("MASTER telescoping identity + closing dl=7 + non-vacuous seed sensitivity", okm and tel==0 and sens)
print("="*70)
npass = sum(1 for _,o in RES if o)
print(f'VERIFY_TOW4R: {npass}/{len(RES)} PASS, total {round(time.time()-t00,1)}s')
for nm,o in RES: print(f'  {nm}: {"PASS" if o else "FAIL"}')
