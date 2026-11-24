
import sys, os, time, math, random, itertools
from fractions import Fraction
import mpmath as mp
import sympy as sp
import flint
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD

print("== (1) split-prime check (Kronecker (D/p)) ==")
expected = {13:[3], 29:[5,7,13], 21:[5], 17:[2,13]}
for D, ps in expected.items():
    got = [p for p in [2,3,5,7,11,13] if sp.kronecker_symbol(D,p)==1]
    print(f"  D={D}: split primes <=13 = {got}  (brief says {ps})  match={got==ps}")
    # inert 2?
    print(f"        (D/2)={sp.kronecker_symbol(D,2)} ({'inert' if sp.kronecker_symbol(D,2)==-1 else 'split' if sp.kronecker_symbol(D,2)==1 else 'ramified'})")

print("== (2) generators: norms, fundamentality of eps via Pell ==")
def min_pell(D):
    # minimal solution x^2 - D y^2 = +-4 with x,y>0, x+y sqrt D minimal; via CF
    # use sympy pell via continued fraction on sqrt(D)
    best=None
    # brute force y up to a bound is fine for these D
    for y in range(1, 2000):
        x2 = D*y*y
        for dlt in (4,-4):
            x = sp.integer_nthroot(x2+dlt,2)
            if x[1]:
                val = Fraction(int(x[0]),2)+Fraction(y,2)*sp.sqrt(D)
                if best is None or val<best[0]:
                    best=(val,int(x[0]),y,dlt)
                break
        if best and Fraction(1,2)*sp.log(best[0]) < sp.log(D*y*y): 
            pass
    return best
for D in (13,29,21,17):
    F = W.make_field(D, split=expected[D][:1] if D!=17 else [2], inert=([] if D==17 else [2]))
    eps = F['unit']
    print(f"  D={D}: eps={W.fmt_qd(eps)} N(eps)={eps.norm()}")
    b = min_pell(D)
    print(f"        minimal Pell (x^2-Dy^2={b[3]}): x={b[1]} y={b[2]} -> ({b[1]}{'+' if b[2]>0 else ''}{b[2]}sqrt{D})/2")
    for i,(g,kind,p) in enumerate(zip(F['letters'],F['kinds'],F['prime_of'])):
        if kind in ('ram','split'):
            print(f"        gen[{F['names'][i]}] = {W.fmt_qd(g)}  norm={g.norm()}  (ideal above {p})")

print("== (3) factorisation test: 30 random S-units per D ==")
random.seed(42)
for D in (13,29,21,17):
    split = {13:[3],29:[5],21:[5],17:[2,13]}[D]
    inert = [] if D==17 else [2]
    F = W.make_field(D, split=split, inert=inert)
    ok=0; fails=[]
    for t in range(30):
        exps=[random.randint(-6,6) for _ in F['letters']]
        x=QD(1,0,D)
        for g,e in zip(F['letters'],exps):
            if e: x=x*(g**e)
        sgn=random.choice([1,-1]); x = x if sgn==1 else -x
        f=W.factor_sunit(x,F)
        if f==exps: ok+=1
        else: fails.append((exps,f))
    print(f"  D={D} (letters {F['names']}): {ok}/30 recovered", "FAILS:"+str(fails[:3]) if fails else "")

print("== (4) planted LLL relation test at 950 digits ==")
mp.mp.dps = 1000
vals = [mp.zeta(5), mp.pi**2*mp.log(2)**3, mp.log(mp.mpf('1.5'))**5, mp.pi**4*mp.log(2)]
a_true = [1048576, -129600, 7776, -349920]  # arbitrary-ish
combo = sum(mp.mpf(ai)*v for ai,v in zip(a_true,vals))
vl = vals + [combo]
rows, dt = W.lll_rows(vl, 950)
found=None
for r in rows:
    a=r[:-1]
    if a[:4]==a_true and a[4]==-1: found=a
print(f"  planted relation recovered: {found is not None} (LLL {dt:.1f}s, dim {len(vl)})")

print("== (5) li5 vs mpmath at 340 digits (planner evaluator) ==")
mp.mp.dps=360
for xs in ('0.3','0.6','0.9','0.94','0.97','0.995','-0.7','-0.95','-0.3'):
    x=mp.mpf(xs)
    a=W.li5(x,340); b=mp.polylog(5,x)
    print(f"  li5({xs}): |planner - mpmath| = {mp.nstr(abs(a-b),3)}")

print("== (6) chi cross-check vs ZQ5 hardcoded ==")
sys.path.insert(0,'../ZQ5_zagier_weight5')
import w5lib as W5
for D in (13,29,21,17,5,8):
    same = all(W.chi(D,a)==W5.chi(D,a) for a in range(0, 5*D))
    print(f"  D={D}: chi match = {same}")
