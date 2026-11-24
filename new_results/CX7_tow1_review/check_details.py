"""Small independent certificates and numerical spot checks for CX7."""
import sys
from pathlib import Path
from math import comb, factorial
import sympy as sp
import mpmath as mp
from sympy.functions.combinatorial.numbers import stirling

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'TOW1_tower_law'))
import towlib as tl

v = sp.Symbol('v')
d, r = sp.symbols('d r')
a, b = 1+sp.sqrt(2), sp.sqrt(2)-1
t = (v*v+2*v-1)/(4*v)
dl1 = 1/(v-b)+1/(v+a)-1/v
dl2 = 1/(v-a)+1/(v+b)-1/v
assert sp.simplify(dl1.subs(v,1/v)/v**2+dl2) == 0
assert sp.simplify(sp.diff(t,v).subs(v,sp.I)) == 0
assert sp.simplify(sp.diff(t,v,2).subs(v,sp.I)) == -sp.I/2
L = (2*d-1)/v+(2-2*d)/(v+sp.I)-2*v/(v*v-a*a)-2*v/(v*v-b*b)
assert sp.simplify(L.subs(v,sp.I)+sp.I*(d-2)) == 0
assert sp.simplify(sp.diff(L,v).subs(v,sp.I)+r*sp.diff(t,v,2).subs(v,sp.I)
                   -((3*d-3)/2-sp.I*r/2)) == 0
assert sp.cancel(-sp.I*(d-1)*(2*d-1)/(2*d-2)+3*sp.I/2+sp.I*(d-2)) == 0
print('PASS symbolic: fold identity, t derivatives, both log derivatives, non-adjacent E1', flush=True)

mp.mp.dps = 80
ii, sq = mp.j, mp.sqrt(2)
aa, bb, ts = 1+sq, sq-1, (1+ii)/2
def bh(n):
    return mp.fprod(mp.mpf('0.5')-j for j in range(n))/factorial(n)
def weight(s,m):
    return mp.mpf(int(stirling(m,s,kind=2)))*factorial(s)*2**m/factorial(m)
def polynomial(D,m):
    terms = [(j,s,bh(j)*comb(j,s)*weight(s,m)*weight(j-s,D-m))
             for j in range(1,D+1) for s in range(j+1)]
    return lambda x: mp.fsum(c*(x*x+(1-x)**2)**(D-j)*x**(2*(j-s))*(1-x)**(2*s)
                             for j,s,c in terms if c)
worst = mp.mpf(0)
for D,m in ((2,2),(3,0),(3,2),(3,3),(4,2),(7,4),(11,8)):
    S = polynomial(D,m)
    value = bh(D)*ii**D*(-1)**m*comb(D,m)
    ratio = 2*(D-2*m)-2*D*ii-mp.mpf(4*(D-1)*(D-2*m))/(2*D-3)
    for actual, expected in ((S(ts),value),(mp.diff(S,ts),value*ratio)):
        err = abs(actual-expected)/max(1,abs(expected))
        worst = max(worst,err)
        assert err < mp.mpf('1e-70'), (D,m,err)
    if m < D-m:
        continue
    def A(z):
        T = (z*z+2*z-1)/(4*z)
        return -sq/2*8**D*z**(2*D-1)*S(T)*(z+ii)**(2-2*D)/((z*z-aa*aa)*(z*z-bb*bb))
    double_factorial = sp.factorial2(2*D-3)
    rho0 = (-1)**D*mp.mpf(int(double_factorial))*ii**(D+1)/(2*sq*factorial(D))
    rho0 *= (-1)**m*(1+(m>D-m))*comb(D,m)
    expected = [rho0, -ii*(D-2)*rho0,
                ((mp.mpf(3*D-3)/2-ii*ratio/2-(D-2)**2)/2)*rho0]
    for delta in range(min(2,2*D-4)+1):
        actual = (1+(m>D-m))*mp.diff(A,ii,delta)/factorial(delta)
        err = abs(actual-expected[delta])/max(1,abs(expected[delta]))
        worst = max(worst,err)
        assert err < mp.mpf('1e-70'), (D,m,delta,err)
print('PASS numerical S, S-prime and rho spot checks: relative-scaled max error', mp.nstr(worst,6), flush=True)

for p in (2,3,4):
    for delta in (0,1,2):
        initial = lambda dp: tl.rho_c((p,p),dp)
        extended = lambda dp: sum((tl.ehat_pp(p,dp,s).scale((-1)**s)
                                  for s in range(dp//2+1)),tl.ZERO)
        assert tl.feed_diag_into_child(p,2*p,delta,initial) == tl.feed_diag_into_child(p,2*p,delta,extended)
    for k in (2*p,2*p+1,2*p+2):
        for delta in (0,1,2):
            assert tl.U_e((p,p),k,delta) == tl.u(p,p,k,delta)
            assert tl.U_e((p,p-1),k,delta) == tl.u(p,p-1,k,delta)
print('PASS exact: E3 and first diagonal/adjacent levels p=2,3,4', flush=True)

for k in (2,3,4):
    assert tl.U_e((1,1),k,0) == tl.u(1,1,k,0)
    h = tl.u(1,1,k,0)
    hp = tl.pf_zero()
    hp[1], hp[2] = [tl.ZERO,h], [tl.ZERO,h.conj_i()]
    G = tl.pf_integrate_higher(hp)
    G[3] = -tl.pf_eval_at(G,tl.BQ)
    feed = tl.pf_add(tl.pf_mul_DL(G,tl.DL2),tl.pf_fold(tl.pf_mul_DL(G,tl.DL1)))
    higher = tl.pf_split(feed)[1]
    assert all(x.is_zero() for part in higher[1:3] for x in part)
print('PASS exact: p=1 degree-two diagonal, k=2,3,4; child degree one has no higher part', flush=True)
