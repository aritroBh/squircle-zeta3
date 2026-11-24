
import os, sys, pickle, time, subprocess
from fractions import Fraction as Fr
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, base + "/research/08_theory/TOW2_tower_all_depths")
from fastred import QSI
H = pickle.load(open(base + "/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl",'rb'))
H.pop('_meta', None)
# my own Q(sqrt2,i) division from scratch (not TOW4's qsi_div, not check_data.Q)
def mydiv(x, y):
    # x,y : 4-tuples (a,b,c,d) = a + b r2 + i(c + d r2). returns ((re,im rational)) if result in Q(i) else (re,im,flag)
    def mul(A,B): return ((A[0]*B[0]+2*A[1]*B[1], A[0]*B[1]+A[1]*B[0]))
    # represent as (R, I) pairs of Q(r2)
    X = ((x[0],x[1]),(x[2],x[3])); Y = ((y[0],y[1]),(y[2],y[3]))
    # 1/Y = conj_s(Y)/(Y conj_s(Y)); Y conj_s(Y) in Q(i)
    Yc = ((Y[0][0],-Y[0][1]),(Y[1][0],-Y[1][1]))
    Nre = mul(Y[0],Yc[0]); Nre = (Nre[0]-mul(Y[1],Yc[1])[0], Nre[1]-mul(Y[1],Yc[1])[1])
    Nim = mul(Y[0],Yc[1]); Nim = (Nim[0]+mul(Y[1],Yc[0])[0], Nim[1]+mul(Y[1],Yc[0])[1])
    assert Nre[1]==0 and Nim[1]==0   # in Q(i)
    den = Nre[0]**2 + Nim[0]**2
    assert den != 0
    inv = ((Nre[0]/den, Fr(0)), (-Nim[0]/den, Fr(0)))
    # X * conj_s(Y) * inv
    def cmul(A,B):
        return ((mul(A[0],B[0])[0]-mul(A[1],B[1])[0], mul(A[0],B[0])[1]-mul(A[1],B[1])[1]),
                (mul(A[0],B[1])[0]+mul(A[1],B[0])[0], mul(A[0],B[1])[1]+mul(A[1],B[0])[1]))
    W = cmul(cmul(X, Yc), inv)
    return W  # ((re_r0, re_r1),(im_r0, im_r1))
import sympy as sp
sys.path.insert(0, work + "/tow4_copy")
import tow4lib as T4
d_, m_ = T4.d, T4.m
def ev_pair(expr, dd, pp):
    v = sp.expand(expr.subs({d_: dd, m_: pp}))
    return Fr(sp.re(v)), Fr(sp.im(v))
t0=time.time()
res = {}
for delta in [7,8]:
    aa = T4.arho(delta)
    bad=0; n=0; sqrt2_seen=0
    for (p,q) in H:
        dd = p+q
        if 2*dd-4 < delta: continue
        num = H[(p,q)][delta][dd]; den = H[(p,q)][0][dd]
        W = mydiv(num, den)
        n += 1
        if W[0][1]!=0 or W[1][1]!=0:
            sqrt2_seen += 1; bad += 1; continue
        er, ei = ev_pair(aa, dd, p)
        if W[0][0] != er or W[1][0] != ei: bad += 1
    res[delta] = (n, bad, sqrt2_seen)
    print(f"RHO my-own-division delta={delta}: {n} cells, bad={bad} (sqrt2-part nonzero in {sqrt2_seen}) ({round(time.time()-t0,1)}s)", flush=True)
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} RHO 717-cell checks with referee's own Q(r2,i) division: {res}\n")
PY = base + "/.venv/bin/python"
r = subprocess.run(["nice","-n","10",PY,"-B","check_full78.py"], capture_output=True, text=True, cwd=work+"/tow4_copy", timeout=1500)
print(r.stdout)
print("stderr:", r.stderr[-200:] if r.stderr else "")
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} reran check_full78.py in tow4_copy: {r.stdout.strip().splitlines()[1]}\n")
