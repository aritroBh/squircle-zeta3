"""Independent initial-data derivation, complete non-adjacent identities, boundary data."""
from pathlib import Path
from math import comb, factorial
from fractions import Fraction as F
import sys
import sympy as sp
import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'TOW2_tower_all_depths'))
import tow2lib as T
import towlib as tl
from genintegrand import S_poly

d,m,v = T.d,T.m,sp.Symbol('v')
i = sp.I
a,b = 1+sp.sqrt(2),sp.sqrt(2)-1
t = v/4+sp.Rational(1,2)-1/(4*v)
assert t.subs(v,i)==(1+i)/2
assert sp.diff(t,v).subs(v,i)==0
for r in range(2,9):
    assert sp.simplify(sp.diff(t,v,r).subs(v,i)-T.tder(r))==0
assert T.check_tau(7)
print('PASS Lemma T derivatives r=1..8 and exact tau table n=0..7',flush=True)

# Direct chain rule: log S(t(v)) contributes R1*t2, R1*t3,
# R1*t4 + 3*(R2-R1^2)*t2^2 at derivative orders 2,3,4.
alg = (2*d-1)/v+(2-2*d)/(v+i)-2*v/(v*v-a*a)-2*v/(v*v-b*b)
ell = {r: sp.simplify(sp.diff(alg,v,r-1).subs(v,i)) for r in range(1,5)}
ell[2] += T.R[1]*sp.diff(t,v,2).subs(v,i)
ell[3] += T.R[1]*sp.diff(t,v,3).subs(v,i)
ell[4] += T.R[1]*sp.diff(t,v,4).subs(v,i)+3*(T.R[2]-T.R[1]**2)*sp.diff(t,v,2).subs(v,i)**2
L = {r: sp.cancel(ell[r]/factorial(r)) for r in ell}
derived = {0:sp.Integer(1),1:L[1],2:L[2]+L[1]**2/2,
           3:L[3]+L[1]*L[2]+L[1]**3/6,
           4:L[4]+L[1]*L[3]+L[2]**2/2+L[1]**2*L[2]/2+L[1]**4/24}
for delta in range(5):
    assert sp.cancel(derived[delta]-T.arho_simple(delta))==0,delta
    assert sp.simplify(T.L(delta)-L[delta])==0 if delta else True
    if delta: print(f'log(A) derivative/{delta}! =',sp.simplify(L[delta]),flush=True)
print('PASS independently derived rho(delta)/rho(0), delta=0..4 equals the submitted cache',flush=True)

def e(delta,s,dd=d,mm=m):
    dd,mm=sp.sympify(dd),sp.sympify(mm)
    if s>delta//2: return sp.Integer(0)
    if s==0: return T.arho_simple(delta).subs({d:dd,m:mm},simultaneous=True)
    if delta==2: return -1/(2*dd-3)
    if delta==3: return i*(dd-3)/(2*dd-3)
    if s==2: return sp.Rational(2)/((2*dd-3)*(2*dd-5))
    return (dd**3-9*dd**2+sp.Rational(113,4)*dd-i*dd/2+i*mm-29)/((2*dd-3)*(2*dd-5))

for delta in range(5):
    for s in range(delta//2+1):
        rhs = sum(((m+1)*T.tau1f(delta-dp)*e(dp,s,d+1,m+1)*T.rA0()
                   +(d-m+1)*T.tau2f(delta-dp)*e(dp,s,d+1,m)*T.rB0())/(2*d-dp-1)
                  for dp in range(delta+1))
        assert sp.cancel(e(delta,s)+e(delta,s+1)-rhs)==0,(delta,s)
        print(f'PASS independent symbolic non-adjacent N({delta},{s}) in d,m',flush=True)

# Small d and endpoint m checks use the actual finite S polynomial, including d<l.
ts=(1+i)/2
for dd in range(1,5):
    for mm in range(dd+1):
        P=S_poly(dd,mm)
        values=[sp.expand(sum(sp.Rational(c)*sp.ff(n,l)*ts**(n-l)
                             for n,c in enumerate(P) if n>=l)) for l in range(5)]
        for l in range(1,5):
            assert sp.simplify(values[l]/values[0]-T.R[l].subs({d:dd,m:mm}))==0,(dd,mm,l)
print('PASS exact S derivative ratios l=1..4 for d=1..4, all m=0..d',flush=True)

# Horner evaluation of the finite polynomial, independently differentiated at 80 digits.
mp.mp.dps=80
ii,sq=mp.j,mp.sqrt(2)
aa,bb,tt=1+sq,sq-1,(1+ii)/2
def num(x): return mp.mpf(str(x.p))/mp.mpf(str(x.q))
def numeric(expr): return mp.mpc(str(sp.N(sp.re(expr),85)),str(sp.N(sp.im(expr),85)))
worst=mp.mpf(0)
samples=[(4,0),(4,2),(4,3),(4,4),(5,3),(7,4),(11,8)]
for dd,mm in samples:
    coeff=[num(sp.Rational(c)) for c in S_poly(dd,mm)]
    def S(z): return mp.polyval(list(reversed(coeff)),z)
    sv=S(tt)
    for l in range(1,5):
        expected=numeric(T.R[l].subs({d:dd,m:mm}))
        err=abs(mp.diff(S,tt,l)/sv-expected)/max(1,abs(expected))
        worst=max(worst,err)
        assert err<mp.mpf('1e-65'),(dd,mm,l,err)
    if mm<dd-mm: continue
    def A(z):
        tv=z/4+mp.mpf('0.5')-1/(4*z)
        return -sq/2*8**dd*z**(2*dd-1)*S(tv)*(z+ii)**(2-2*dd)/((z*z-aa*aa)*(z*z-bb*bb))
    for delta in (3,4):
        actual=(1+(mm>dd-mm))*mp.diff(A,ii,delta)/factorial(delta)
        exact=tl.u(mm,dd-mm,dd,delta).sympy()
        expected=numeric(exact)
        err=abs(actual-expected)/max(1,abs(expected))
        worst=max(worst,err)
        assert err<mp.mpf('1e-65'),(dd,mm,delta,err)
        assert sp.simplify(exact-tl.u(mm,dd-mm,dd,0).sympy()*T.arho_simple(delta).subs({d:dd,m:mm}))==0
print('PASS numerical S derivatives through 4 and rho(3),rho(4) vs exact initial towers; max scaled error',mp.nstr(worst,8),flush=True)

def eqsi(expr):
    z=sp.expand(expr)
    return tl.QSI(F(sp.Rational(sp.re(z))),0,F(sp.Rational(sp.im(z))),0)
def off(key,delta,s):
    return tl.rho(key,0)*eqsi(e(delta,s,sum(key),key[0]))
def hat(p,delta,s):
    return sum((off((p+1,p),dp,s)*tl.tau1(delta-dp).scale(F(p+1,4*p-dp-1))
                for dp in range(delta+1)),tl.ZERO)
def candidate(key,k,delta):
    p,q=key
    if p==q:
        if k==2*p:return tl.rho(key,delta)
        return sum((hat(p,delta,s).scale(comb(k-2*p-1,s)) for s in range(delta//2+1)),tl.ZERO)
    return sum((off(key,delta,s).scale(comb(k-p-q,s)) for s in range(delta//2+1)),tl.ZERO)
count=0
for p in (1,2,3,4):
    for key in ((p,p),(p,p-1)):
        if sum(key)<2:continue
        for k in range(max(sum(key),2*p-1),2*p+4):
            for delta in range(min(4,2*sum(key)-4)+1):
                assert candidate(key,k,delta)==tl.u(*key,k,delta),(key,k,delta)
                count+=1
    if p>=2:
        for delta in range(min(4,4*p-6)+1):
            initial=lambda dp: tl.rho((p,p),dp)
            extended=lambda dp:sum((hat(p,dp,s).scale((-1)**s) for s in range(dp//2+1)),tl.ZERO)
            assert tl.feed_diag_into_child(p,2*p,delta,initial)==tl.feed_diag_into_child(p,2*p,delta,extended),(p,delta)
print(f'PASS {count} exact diagonal/adjacent boundary cells p=1..4; E3 feeds at every admissible depth <=4 for p=2,3,4',flush=True)

# Ensure the independently certified N_l polynomials are exactly the library's candidates.
for name,stop in [('prove_S12.py','for l in (1, 2):'),('prove_S34.py','import sys, time')]:
    ns={}
    exec((HERE.parent/'TOW2_tower_all_depths'/name).read_text().split(stop)[0],ns)
    for l,numerator in ns['N'].items():
        assert sp.cancel(numerator-T.N[l])==0,(name,l)
print('PASS fitted N1..N4 in submitted certificates equal the library polynomials',flush=True)

# Independent coefficient derivation: set h=t-t*, q2=2*i*h+2*h^2,
# j=d-n and s=m-r. Cancel C(1/2,d-n)/C(1/2,d) / d_falling_n
# to (-2)^n/Pi_n(d), before using any endpoint values.
x=sp.Symbol('x')
tails={0:sp.Integer(1),1:x*(x-1)/2,2:x*(x-1)*(x-2)*(3*x-5)/24,
       3:x*(x-1)*(x-2)**2*(x-3)**2/48,
       4:x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760}
zs=1-ts
for l in (3,4):
    result=0
    for n in range(l+1):
        for r in range(n+1):
            base=(-2)**n/T.Pi(n)*tails[r].subs(x,m)*tails[n-r].subs(x,d-m)*ts**(2*(r-n))*zs**(-2*r)
            aa,bb=2*(d-m-n+r),2*(m-r)
            for j in range(min(n,l-n)+1):
                rem=l-n-j
                bc=sum(sp.ff(aa,c)/factorial(c)/ts**c*sp.ff(bb,rem-c)/factorial(rem-c)*(-1/zs)**(rem-c)
                       for c in range(rem+1))
                result+=base*2**n*comb(n,j)*i**(n-j)*bc
    assert sp.cancel(sp.expand(factorial(l)*result)-T.R[l])==0,l
    print(f'PASS independent h-coefficient derivation S^{l}/S for symbolic d,m',flush=True)

p=sp.Symbol('p',integer=True,positive=True)
def normalized_hat(delta,s):
    return (p+1)*sum(T.tau1f(delta-dp)*e(dp,s,2*p+1,p+1)*2*i*(4*p-1)/(p+1)/(4*p-dp-1)
                     for dp in range(delta+1))
claims={3:-(2*p-3)*(32*p**3-96*p**2+130*p-69)/(12*(4*p-3)),
        4:i*(512*p**6-3840*p**5+12992*p**4-25600*p**3+30722*p**2-20789*p+6054)/(48*(4*p-5)*(4*p-3))}
for delta in (3,4):
    discrepancy=T.arho_simple(delta).subs({d:2*p,m:p})-sum((-1)**s*normalized_hat(delta,s) for s in range(delta//2+1))
    assert sp.cancel(discrepancy-claims[delta])==0
    print(f'PASS displayed Delta({delta}) is normalized by rho_pp(0); diagonal leading Newton ratio =',sp.factor(normalized_hat(delta,delta//2)),flush=True)
