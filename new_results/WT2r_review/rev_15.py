
import sys, pickle
sys.path.insert(0,'.')
import sympy as sp
from fastred import *
from fastred import v
from fractions import Fraction as Fr
from math import comb, factorial
record = pickle.load(open('record_W2b.pkl','rb'))
S2=sp.sqrt(2); I=sp.I
tau1 = QSI(-1,0,-1,0)/QSI(2); tau2 = QSI(-1,0,1,0)/QSI(2)
IQ=QSI(0,0,1,0); SQ=QSI(0,1,0,0)
def Kf(d):
    odd=1
    for m in range(1,2*d-2,2): odd*=m
    return (IQ**(d+1))*QSI(Fr((-1)**d*odd,2*factorial(d)))/SQ
def rho_cf(p,q):
    d=p+q
    return Kf(d)*((-ONE)**p)*QSI(comb(d,p)*(1+(1 if p>q else 0)))
def h(k,key):
    i,j=key; d=i+j
    rec=record.get((k,key))
    if rec is None: return None
    H=rec['H'][1]; top=2*d-2
    return H[top-1] if len(H)>=top else ZERO
def h_hat(p): return rho_cf(p+1,p)*QSI(p+1)*tau1/QSI(4*p-1)

# 1. rho vs initial data (reducer), d<=12
ok = all((record[(d,(p,q))]['H'][1][2*d-3] if len(record[(d,(p,q))]['H'][1])>=2*d-2 else ZERO)==rho_cf(p,q)
         for d in range(2,13) for p in range(d,(d-1)//2,-1) for q in [d-p])
print("1. rho_pq closed form == reducer initial data, d=2..12:", ok)

# 2. Lemma D stabilization vs all data
okD=True
for k in range(2,13):
    for i in range(1,k+1):
        for j in range(0,i+1):
            d=i+j
            if d<2 or d>k: continue
            hv=h(k,(i,j))
            if hv is None: continue
            pred = rho_cf(i,j) if i>j else (rho_cf(i,i) if k==2*i else h_hat(i))
            okD &= (hv==pred)
print("2. Lemma D (h=rho off-diag; h_pp=rho_pp at k=2p, hat_p for k>2p), all data k<=12:", okD)

# 3. HigherPart of (2,0),(1,1) has length exactly <=2 for all levels (top IS everything)
lens = {(k,key): (len(record[(k,key)]['H'][1]), len(record[(k,key)]['H'][2])) for k in range(2,13) for key in [(2,0),(1,1)]}
print("3. max +-i order of HigherPart((2,0)),((1,1)) over k=2..12:", max(x[0] for x in lens.values()), "(bound 2d-2=2)")

# 4. g20, g11 from rho's: Gm = (1-eval_b) Int[h u^-2 + conj]
def gm_from_h(hv):
    # HigherPart = hv/(v-i)^2 + conj(hv)/(v+i)^2 ; Int = -hv/(v-i) - conj/(v+i) ; subtract value at b
    hs = hv.sympy(); hc = hv.conj_i().sympy()
    b = S2-1
    G = -hs/(v-I) - hc/(v+I)
    Gb = G.subs(v, b)
    return sp.radsimp(sp.simplify(G - Gb))
g20 = (1+S2)/4 - S2/(2*(v**2+1))
g11 = -S2*(v**2+2*v-1)/(4*(v**2+1))
d20 = sp.simplify(gm_from_h(rho_cf(2,0)) - g20)
d11a = sp.simplify(gm_from_h(rho_cf(1,1)) - (-g20))
d11b = sp.simplify(gm_from_h(h_hat(1)) - g11)
print("4. Gm from rho_20 == g20:", d20==0, "; Gm from rho_11 == -g20:", d11a==0, "; Gm from hat_1 == g11:", d11b==0)
print("   rho_20 =", rho_cf(2,0).sympy(), " rho_11 =", rho_cf(1,1).sympy(), " hat_1 =", h_hat(1).sympy())

# 5. Gm data check
g20_pf = pf_from_sympy(g20); g11_pf = pf_from_sympy(g11); ng20_pf = pf_neg(g20_pf)
def pf_nz(P):
    return bool(P[0]) or any(not e.is_zero() for e in P[1]) or any(not e.is_zero() for e in P[2]) or not P[3].is_zero()
okG = all(not pf_nz(pf_sub(record[(k,(2,0))]['Gm'], g20_pf)) for k in range(2,13))
okG11 = all(not pf_nz(pf_sub(record[(k,(1,1))]['Gm'], ng20_pf if k==2 else g11_pf)) for k in range(2,13))
print("5. Gm_20==g20 (k=2..12):", okG, " Gm_11==-g20 (k=2), ==g11 (k>=3):", okG11)

# 6. SimplePart(W_10)==Rw2 and eps-identity
Rw2   = 1/(2*v) + 1/(2*(v-(1+S2))) - 1/(v+(1+S2)) - 1/(2*(v+(S2-1)))
Rw2_pf = pf_from_sympy(Rw2)
okW = all(not pf_nz(pf_sub(pf_split(record[(k,(1,0))]['W'])[0], Rw2_pf)) for k in range(2,13))
print("6. SimplePart(W_10^(k)) == Rw2, k=2..12:", okW)
eps = sp.Symbol('eps')
a=1+S2; b=S2-1
R10   = (1-S2/2)/v - 2/(v+a)
DL1s = 1/(v-b) + 1/(v+a) - 1/v
DL2s = 1/(v-a) + 1/(v+b) - 1/v
def FF(R): return sp.together(R.subs(v, 1/v)/v**2)
Gm20 = eps**2/(1-eps)*g20
Gm11 = -eps**2*g20 + eps**3/(1-eps)*g11
Wcas = -2*Gm20*DL1s - Gm11*DL2s - FF(sp.together(Gm11*DL1s))
num = sp.expand(sp.together((1-eps)*Wcas))
c2 = sp.together(num.coeff(eps,2)); c3 = sp.together(num.coeff(eps,3))
rem = sp.together(num - c2*eps**2 - c3*eps**3)
ok6 = (rem==0) and sp.simplify(sp.radsimp(sp.together(c2-Rw2)))==0 and sp.simplify(sp.radsimp(sp.together(c3)))==0
print("7. (1-eps)(W_10-eps R_10) = eps^2 Rw2 exactly (c2=Rw2, c3=0, no other terms):", ok6)
ap = sp.apart(sp.together(Wcas), v, extension=[S2, I])
ipoles = [tt for tt in sp.Add.make_args(sp.expand(ap)) if sp.denom(sp.together(tt)).has(I)]
print("8. poles at +-i surviving in W_10(eps):", ipoles if ipoles else "NONE")
