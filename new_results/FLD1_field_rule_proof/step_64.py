WD='<project>/research/08_theory/FLD1_field_rule_proof'
ver = r'''"""verify_fld1.py -- self-contained verification of the Field-rule proof (FLD1), k=2,3 reproduction + k=4 checks.

Pipeline (all exact unless stated):
  [1] symbolic-n mu_{k,n} (Lemma A) -> Z-sum canonical form (k<=4)   [step_22/23 machinery]
  [2] CONSTRUCT the HPL forms (Lemma C): bijection recursion + symmetrize + correction, then PROVE
      for all n:  even Z-sum identity [c^(2n)]F == mu  AND  odd [c^(2n+1)]F == 0  AND  low orders c^-3..c^3.
  [3] k=2,3 reproduction: constructed forms == fitted B2 (thm:M2alpha) and B3 (Taylor to c^40).
  [4] k=4 checks: exact [c^(2n)]F == mu_{4,n} for n<=100 (all zeta-monomials), odd/non-positive = 0;
      and NUMERICAL agreement with the corner integral (corner_series) at 3 angles to >=40 digits.
Run:  python verify_fld1.py      (expected < 5 min)
"""
import sys, time, json, re, pickle
from fractions import Fraction as Fr
sys.path.insert(0, '.')
import sympy as sp
t0 = time.time()

def exectop(fn, argv1=None):
    old = sys.argv
    if argv1 is not None: sys.argv = ['x', argv1]
    src = open(fn).read(); ns = {}
    exec(src[:src.index("if __name__")], ns)
    sys.argv = old
    return ns

PASS = []
def check(name, ok):
    PASS.append((name, bool(ok)))
    print(f'  [{"PASS" if ok else "FAIL"}] {name}   (t={time.time()-t0:.1f}s)', flush=True)

print("="*70)
print("[1] Lemma A: symbolic mu_{k,n} and Z-sum canonical form (k<=4)")
ns22 = exectop('step_22_musymbolic.py', '4')   # KMAX=4
import zsum
from zsum import K, nn, ONE, ZERO, add_into, lin
mu_formula = ns22['mu_formula']; KMAX = ns22['KMAX']
def harm_n_minus_1(r): return {(r,): K(2**(r-1)), (-r,): K(2**(r-1)), (): -ONE/nn**r}
def harm_2n_minus_1(r): return {(r,): ONE, (): -ONE/(2*nn)**r}
def to_zsum(F):
    parts = {}
    for mon, cf in F.terms():
        E = {(): ONE}
        for r in range(1, KMAX+1):
            for _ in range(mon[r-1]): E = zsum.mul(E, harm_n_minus_1(r))
            for _ in range(mon[KMAX+r-1]): E = zsum.mul(E, harm_2n_minus_1(r))
        zkey = tuple(list(mon[2*KMAX:])+[0]*(4-len(mon[2*KMAX:])))[:4]
        parts.setdefault(zkey, {})
        for w, c in E.items(): add_into(parts[zkey], w, cf*c)
    return parts
Zin = {k: to_zsum(mu_formula(k)) for k in range(KMAX+1)}
# verify mu symbolic against fast numeric table (n<=60) -- reuse step_21's mu_exact json
MU = json.load(open('mu_exact_K5_N100.json'))
def Hf(r,m): return sum((Fr(1,i**r) for i in range(1,m+1)),Fr(0))
def evalK(c,nv): return Fr(int(c.numer(nv)),1)/Fr(int(c.denom(nv)),1)
okA = True
for k in range(KMAX+1):
    F = mu_formula(k)
    for nv in range(1, 61):
        acc={}
        for mon, cf in F.terms():
            cval=evalK(cf,nv); hv=Fr(1)
            for r in range(1,KMAX+1):
                hv*=Hf(r,nv-1)**mon[r-1]; hv*=Hf(r,2*nv-1)**mon[KMAX+r-1]
            ze=tuple(list(mon[2*KMAX:])+[0]*(4-len(mon[2*KMAX:])))[:4]
            zk=','.join(map(str,ze)); acc[zk]=acc.get(zk,Fr(0))+cval*hv
        ref=MU[str(k)][str(nv)]
        for kk in set(list(acc)+list(ref)):
            if Fr(acc.get(kk,'0'))!=Fr(ref.get(kk,'0')): okA=False
check('Lemma A: symbolic mu_{k,n} == exact numeric table (k<=4, n<=60)', okA)

print("="*70)
print("[2] Lemma C: construct HPL forms and PROVE all-n identities (k<=4)")
ns25 = exectop('step_25_construct.py', '3')
construct=ns25['construct']; symmetrise=ns25['symmetrise']; correct=ns25['correct']
even_form=ns25['even_form']; odd_form=ns25['odd_form']; toK=ns25['toK']
from hpl import coeffs_exact
def low_order(terms, k, zkey):
    zks=','.join(map(str,zkey))
    for N in range(-3,4):
        tot=Fr(0)
        for (e,w,q) in terms:
            if len(w)==0:
                if N==e: tot+=q
                continue
            idx=N-e
            if idx>=1: tot+=q*coeffs_exact(w,idx+2)[idx]
        target=Fr(0)
        if N>=2 and N%2==0: target=Fr(MU[str(k)][str(N//2)].get(zks,'0'))
        if tot!=target: return False
    return True
store={}
okC=True
for k in range(KMAX+1):
    for zkey, E in sorted(Zin[k].items()):
        Elist=[(w, c.as_expr()) for w,c in E.items()]
        terms=construct(Elist); Fsym=symmetrise(terms); Fsym,bad=correct(Fsym,k)
        store[(k,zkey)]=Fsym
        ev=even_form(Fsym); od=odd_form(Fsym)
        muE={tuple(w): toK(c) for w,c in Elist}
        diff=lin((1,ev),(-1,muE))
        lo=low_order(Fsym,k,zkey)
        if not(len(diff)==0 and len(od)==0 and lo): okC=False; print('   problem at',k,zkey,len(diff),len(od),lo)
check('Lemma C: all-n even identity [c^2n]F==mu AND odd==0 AND low orders, for k<=4 all zeta-monomials', okC)

print("="*70)
print("[3] k=2,3 reproduction: constructed == fitted B2 (thm:M2alpha) and B3")
G1='<project>/research/10_frontier/G1_zeta3_any_angle/'
def load_fit(path):
    out=[]
    for nm,v in json.load(open(path)):
        v=Fr(v)
        if nm=='1': out.append((0,(),v)); continue
        if nm=='c': out.append((1,(),v)); continue
        m=re.match(r'c\^(-?\d+) H\[(.*)\]', nm); out.append((int(m.group(1)),tuple(int(t) for t in m.group(2).split(',')),v))
    return out
def ser(terms,lo,hi):
    out={}
    for (e,w,q) in terms:
        if len(w)==0:
            if lo<=e<=hi: out[e]=out.get(e,Fr(0))+q
            continue
        ce=coeffs_exact(w,hi-e+2)
        for j in range(lo,hi+1):
            idx=j-e
            if 1<=idx<len(ce): out[j]=out.get(j,Fr(0))+q*ce[idx]
    return out
FIT={(2,(0,0,0,0)):'hpl_fit_k2_(0,0,0).json',(3,(0,0,0,0)):'hpl_fit_k3_(0,0,0).json',
     (3,(1,0,0,0)):'hpl_fit_k3_(1,0,0).json',(3,(0,1,0,0)):'hpl_fit_k3_(0,1,0).json'}
ok3=True
for key,fn in FIT.items():
    fit=load_fit(G1+fn); mine=store.get(key,[])
    sf=ser(fit,-3,40); sm=ser(mine,-3,40)
    diff={j:sf.get(j,Fr(0))-sm.get(j,Fr(0)) for j in set(sf)|set(sm)}
    diff={j:v for j,v in diff.items() if v!=0}
    if diff: ok3=False; print('   mismatch',key,list(diff.items())[:3])
check('k=2 construction == thm:M2alpha (B2); k=3 construction == B3 (Taylor to c^40)', ok3)

print("="*70)
print("[4] k=4 checks: exact mu (n<=100) + corner_series numerical")
# exact [c^N] for N=-4..200 vs mu, all zeta-monomials
ok4=True
for (k,zkey),terms in store.items():
    if k!=4: continue
    zks=','.join(map(str,zkey))
    for N in range(-4,201):
        tot=Fr(0)
        for (e,w,q) in terms:
            if len(w)==0:
                if N==e: tot+=q
                continue
            idx=N-e
            if idx>=1: tot+=q*coeffs_exact(w,idx+2)[idx]
        target=Fr(0)
        if N>=2 and N%2==0: target=Fr(MU[str(k)][str(N//2)].get(zks,'0'))
        if tot!=target: ok4=False; print('   exact mismatch',k,zkey,N)
check('k=4 exact [c^2n]F==mu_{4,n} for n<=100, odd/non-positive==0', ok4)

# numerical vs corner_series
import mpmath as mp
mp.mp.dps = 60
sys.path.insert(0,'<project>/research/10_frontier/G1_zeta3_any_angle/')
import corner_series
from hpl import coeffs_mp
def zeta_of(zkey):
    v=mp.mpf(1)
    for i,e in enumerate(zkey):
        if e: v*=mp.zeta(i+2)**e
    return v
N=500
words=set()
for (k,zk),terms in store.items():
    for (e,w,q) in terms:
        if len(w)>0: words.add(w)
wc={w: coeffs_mp(w,N) for w in words}
def Mk_value(k,c):
    tot=mp.mpf(0)
    for (kk,zk),terms in store.items():
        if kk!=k: continue
        s=mp.mpf(0)
        for (e,w,q) in terms:
            qq=mp.mpf(q.numerator)/q.denominator
            if len(w)==0: s+=qq*c**e; continue
            ce=wc[w]; sv=mp.mpf(0); xp=c
            for m in range(1,N+1): sv+=ce[m]*xp; xp*=c
            s+=qq*c**e*sv
        tot+=zeta_of(zk)*s
    return tot
oknum=True
for tag,a in [('2pi/3',2*mp.pi/3),('pi/2',mp.pi/2),('3pi/4',3*mp.pi/4)]:
    c=mp.cos(a/2); Mc=corner_series.corner_coeffs(a,4)
    d=abs(Mk_value(4,c)-Mc[4])
    digits = -mp.log10(d+mp.mpf('1e-999'))
    print(f'   {tag} (c={mp.nstr(c,6)}): |HPL form - corner_series M4| = {mp.nstr(d,3)}  (~{mp.nstr(digits,4)} digits)')
    if d > mp.mpf('1e-40'): oknum=False
check('k=4 HPL form == corner_series M4 at 3 angles to >=40 digits', oknum)

print("="*70)
n_pass = sum(1 for _,o in PASS if o)
print(f'RESULT: {n_pass}/{len(PASS)} checks passed in {time.time()-t0:.1f}s')
for nm,o in PASS: print(('  PASS ' if o else '  FAIL ')+nm)
sys.exit(0 if n_pass==len(PASS) else 1)
'''
open(WD+'/verify_fld1.py','w').write(ver)
import subprocess
VENV='<project>/.venv/bin/python'
import time
t=time.time()
r = subprocess.run(['nice',VENV,'verify_fld1.py'],capture_output=True,text=True,cwd=WD)
print("TOTAL time", round(time.time()-t,1),"s")
print(r.stdout); print(r.stderr[-3000:] if r.returncode else '')