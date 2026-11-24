"""step31_k6.py -- k=6: symbolic mu -> Z-sum -> verify vs table, properness, weight budget."""
import sys, json, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import sympy as sp
import zsum
from zsum import K, nn, ONE, ZERO, add_into
sys.argv=['x','6']
src=open("step_22_musymbolic.py").read(); ns={}
exec(src[:src.index("if __name__")], ns)
KMAX=ns['KMAX']; mu_formula=ns['mu_formula']
MU=json.load(open('mu_exact_K6_N60.json'))
t0=time.time()
print("computing mu_formula(6)...", flush=True)
F6=mu_formula(6)
print("  done", round(time.time()-t0,1),"s,", len(F6.terms()),"monomials", flush=True)
# verify vs table
def Hf(r,m): return sum((Fr(1,i**r) for i in range(1,m+1)),Fr(0))
def evalK(c,nv): return Fr(int(c.numer(nv)),1)/Fr(int(c.denom(nv)),1)
ok=True
for nv in range(1,61):
    acc={}
    for mon,cf in F6.terms():
        cval=evalK(cf,nv); hv=Fr(1)
        for r in range(1,KMAX+1):
            hv*=Hf(r,nv-1)**mon[r-1]; hv*=Hf(r,2*nv-1)**mon[KMAX+r-1]
        ze=tuple(list(mon[2*KMAX:])+[0]*(5-len(mon[2*KMAX:])))[:5]
        zk=','.join(map(str,ze)); acc[zk]=acc.get(zk,Fr(0))+cval*hv
    ref=MU['6'][str(nv)]
    for kk in set(list(acc)+list(ref)):
        if Fr(acc.get(kk,'0'))!=Fr(ref.get(kk,'0')): ok=False; print('MISMATCH',nv,kk)
print("mu_formula(6) == numeric table (n<=60):", ok, flush=True)
# to_zsum
def harm_n_minus_1(r): return {(r,): K(2**(r-1)), (-r,): K(2**(r-1)), (): -ONE/nn**r}
def harm_2n_minus_1(r): return {(r,): ONE, (): -ONE/(2*nn)**r}
def to_zsum(F):
    parts={}
    for mon,cf in F.terms():
        E={(): ONE}
        for r in range(1,KMAX+1):
            for _ in range(mon[r-1]): E=zsum.mul(E,harm_n_minus_1(r))
            for _ in range(mon[KMAX+r-1]): E=zsum.mul(E,harm_2n_minus_1(r))
        zkey=tuple(list(mon[2*KMAX:])+[0]*(5-len(mon[2*KMAX:])))[:5]
        parts.setdefault(zkey,{})
        for w,c in E.items(): add_into(parts[zkey],w,cf*c)
    return parts
t=time.time()
Z=to_zsum(F6)
print("to_zsum:", round(time.time()-t,1),"s; words per zeta-mono:", {zk:len(E) for zk,E in Z.items()}, flush=True)
pickle.dump({zk:[(w,c.as_expr()) for w,c in E.items()] for zk,E in Z.items()}, open('mu_zsum_K6.pkl','wb'))
# properness + weight budget
n=sp.Symbol('n')
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))
def analyze(zkey, terms):
    maxdeg=-999; poles=set(); worst=0
    for w,cs in terms:
        e=sp.together(sp.sympify(cs)); num,den=sp.fraction(e)
        dd=sp.Poly(sp.expand(den),n); dn=sp.Poly(sp.expand(num),n)
        maxdeg=max(maxdeg, dn.degree()-dd.degree())
        if dd.degree()>0:
            rts=sp.roots(sp.expand(den),n)
            for rt,mult in rts.items():
                poles.add(sp.nsimplify(rt))
            worst=max(worst, sum(abs(x) for x in w)+max(rts.values()))
        else:
            worst=max(worst, sum(abs(x) for x in w))
    return maxdeg, poles, worst
print("k=6 analysis:")
allok=True
for zkey,terms in sorted(Z.items()):
    md,poles,worst=analyze(zkey,terms)
    bound=6+1-wt(zkey)
    okb = worst<=bound
    allok = allok and okb and (md<=-2)
    print(f"  z={zkey} (wt {wt(zkey)}): deg={md}, max(|w|+po)={worst} (bound {bound}) ok={okb and md<=-2}; poles={sorted(poles,key=sp.N)}")
print("k=6 proper (deg<=-2) AND weight budget (<=k+1-wt(z)):", allok)
