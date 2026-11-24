import sys
sys.path.insert(0,'.')
import sympy as sp
from fastred import *
from fastred import v
import genintegrand, buildintegrand
from genintegrand import S_poly, w, binom_half, p_comp_1mt
from buildintegrand import sqrt_part_PF, Stilde
from math import comb
S2=sp.sqrt(2); I=sp.I

# The analytic formula for the unfolded sqrt part at key (m, d-m):
# R~(v) = (1/(2 sqrt2)) sum_{j=1}^d binom(1/2,j) z^{2-2j} U_j(t) / (v t (1-t))
# with z=(v^2+1)/(2 sqrt2 v), t=(v^2+2v-1)/(4v), U_j(t)=sum_s binom(j,s) w(s,m) w(j-s,d-m) (1-t)^{2s} t^{2(j-s)}
t_ = (v**2+2*v-1)/(4*v)
z_ = (v**2+1)/(2*S2*v)
def Uj_poly(j,m,d):
    # as sympy rational in t
    tt = sp.Symbol('T')
    out = 0
    for s in range(0,j+1):
        w1=w(s,m); w2=w(j-s,d-m)
        if w1==0 or w2==0: continue
        out += comb(j,s)*sp.Rational(w1.numerator,w1.denominator)*sp.Rational(w2.numerator,w2.denominator)*(1-tt)**(2*s)*tt**(2*(j-s))
    return sp.expand(out)

def Rtilde(d,m):
    tt=sp.Symbol('T')
    ssum=0
    for j in range(1,d+1):
        bh=binom_half(j)
        U=Uj_poly(j,m,d)
        ssum += sp.Rational(bh.numerator,bh.denominator)*sp.expand(U) * z_**(2-2*j)
    ssum = sp.expand(ssum).subs(tt, t_)
    return sp.together(ssum/(2*S2*v*t_*(1-t_)))

# compare with buildintegrand sqrt_part_PF(k, S_m)
bad=[]
for d in range(2,7):
    for m in range(0,d+1):
        S = S_poly(d,m)
        PF = sqrt_part_PF(d,S)
        a = sp.together(pf_to_sympy(PF))
        b_ = Rtilde(d,m)
        if sp.simplify(sp.together(a-b_))!=0:
            bad.append((d,m))
print("analytic R~ formula == buildintegrand sqrt_part (d<=6, all m):", "ALL MATCH" if not bad else bad)

# folded piece formula check: buildintegrand folds at t-level: sqrt_part with S_q(1-t).
# its v-level meaning: pf_fold of the (q,p) piece. Verify pf_fold(sqrt_part(S_q)) == sqrt_part(S_q(1-t)) for small d,q
bad=[]
for d in range(2,7):
    for q in range(0,(d)//2):
        p=d-q
        PFq = sqrt_part_PF(d, S_poly(d,q))
        folded = pf_fold(PFq)
        PFbuilt = sqrt_part_PF(d, p_comp_1mt(S_poly(d,q)))
        dff = pf_sub(folded, PFbuilt)
        nz = (dff[0] and any(not c.is_zero() for c in dff[0].values())) or any(not e.is_zero() for e in dff[1]) or any(not e.is_zero() for e in dff[2]) or not dff[3].is_zero()
        if nz: bad.append((d,q))
print("pf_fold(sqrt part S_q) == sqrt_part with S_q(1-t) (d<=6):", "ALL MATCH" if not bad else bad)

# symbolic lemmas: z'(i)=sqrt2/2, C_d^{(m)} = 2^d (-1)^m binom(d,m), U_d^{(m)}(tau)=i^d(-1)^m binom(d,m)
print("z'(i) =", sp.simplify(sp.diff(z_,v).subs(v,I)), " (expect sqrt(2)/2)")
bad=[]
for d in range(2,9):
    for m in range(0,d+1):
        U=Uj_poly(d,m,d)
        val = sp.expand(U).subs(sp.Symbol('T'), (1+I)/2)
        pred = I**d * (-1)**m * comb(d,m)
        if sp.simplify(val-pred)!=0: bad.append((d,m,val,pred))
print("U_d^{(m)}((1+i)/2) == i^d (-1)^m binom(d,m)  (d<=8):", "ALL HOLD" if not bad else bad[:4])
