"""RST1R step_06: R8 monotonicity — analytic closure for ALL k>=99 + certified grid.
Pieces normalized by t(k)=(k+1)^(-1/2) e^(-sqrt(2 pi (k+1))).
Part A: certified P(k+1)<=P(k) on k in [99,3000] (arb balls, ratio upper < 1).
Part B: analytic closures (arb-certified coefficient inequalities):
  pairing trick: -A/sqrt(j) + B/sqrt(N) + 1/(2N) + R(k) < 0 for all k>=99,
  using 1/(2N) <= 0.05/sqrt(N), (2k+1)/(2k(k+1)) <= 0.10051/sqrt(N) (k>=99).
"""
from flint import arb, ctx
ctx.dps = 50
PI = arb.pi(); SQ2 = arb(2).sqrt(); sqpi = PI.sqrt(); sq2pi=(2*PI).sqrt()
a = PI/(2*SQ2)

def FH_gauss(j, q, lam):
    j,q = arb(j),arb(q)
    rL=(j*q/(SQ2*lam)).sqrt(); c0=1+q/SQ2; Sp=(2*SQ2*lam*j*q).sqrt()
    CL=2**arb('1.5')/(3*lam); c2=8*SQ2*j*q/(4*rL+c0)**3
    mid=((rL-c0)/(2*SQ2))**arb('-1.5')*(2*PI/c2).sqrt()
    return sqpi/q*(CL*(-Sp/4).exp()+mid+(SQ2*rL)**arb('-1.5')/lam)*(lam*c0/2).exp()*(-Sp).exp()

def tg(k):
    k=arb(k)
    return (k+1)**arb('-0.5')*(-sq2pi*(k+1).sqrt()).exp()

def P_q2(k):   return 8/PI*(FH_gauss(k,2,a)+2*FH_gauss(k-1,2,a))/tg(k)
def P_mid(k):  return sum(8/PI*(arb(2)**i/arb(i).fac())*FH_gauss(k-i,i+1,a)/tg(k) for i in range(2,13))
def P_mid2(k):
    k=arb(k); Sinf = arb(1).exp()**2 - sum(arb(2)**i/arb(i).fac() for i in range(0,13))
    c013 = 1+14/SQ2
    return 8/PI*2*SQ2/14*(2*(k*(k+1)/PI).sqrt()*(a*c013/2).exp() + 1/a)*(-(14*PI*(k-13)).sqrt()).exp()*Sinf/tg(k)
def P_tail(k, fixed=True):
    k=arb(k); xp=2*(a/(2*SQ2)).exp()
    s1 = (a/(2*SQ2)).exp()*(arb(1).exp()*xp).exp()
    if fixed: s1 = (a/2).exp()*s1
    s2 = (2*arb(1).exp()).exp()
    return 8/PI*2*SQ2/(k/2+1)*(2*(k*(k+1)/PI).sqrt()*s1 + s2/a)*(-k/2).exp()/tg(k)
def P_T2(k):   return 4/PI*FH_gauss(k,1,3*a)/tg(k)
def P_T3J(k):
    k=arb(k); return (2/(PI*k)*arb(3)**(-k) + arb(2)**(-k))/tg(k)

print("== Part A: certified grid P(k+1)<=P(k) on [99,3000] (arb) ==")
pieces = [('q2',P_q2),('mid',P_mid),('mid2',P_mid2),('tail_fixed',P_tail),('T2',P_T2),('T3J',P_T3J)]
for name, f in pieces:
    worst = arb(0); wk=None
    prev = f(99)
    for k in range(99, 3000):
        cur = f(k+1)
        ratio = cur/prev
        if ratio.upper() > worst: worst = ratio.upper(); wk=k
        prev = cur
    print(f"  {name:10s}: max P(k+1)/P(k) = {worst.str(10)} at k={wk}  ({'<1 OK' if worst<1 else 'FAIL'})", flush=True)
    assert worst < 1

print("== Part B: analytic closures, certified ==")
B = sq2pi/2
one99 = arb(1)/arb(99+1).sqrt()
print(f"  k>=99: 1/(2N) <= (1/(2 sqrt(N)))*sqrt(N)/N... using 1/sqrt(N)<=0.1: 1/(2N) <= 0.05/sqrt(N); (2k+1)/(2k(k+1)) <= 0.10051/sqrt(N)")
# q2-left (j=k-1 worst): -(5/4) sq2pi/(2 sqrt j) + sq2pi/(2 sqrt N) + 1/(2N) < 0
# suffices: (5/4)B^2*(k+1) > (B+0.05)^2*(k-1)  [since 5B/(4 sqrt j) >= 5B/(4 sqrt N)]
A = arb(5)/4*B
print(f"  q2-left: A^2={ (A**2).str(8)} > (B+0.05)^2={((B+arb('0.05'))**2).str(8)}; linear-in-k inequality A^2(k+1)>(B+0.05)^2(k-1): coeff {((A**2-(B+arb('0.05'))**2)).str(6)} > 0 => holds all k>=99")
assert A**2 > (B+arb('0.05'))**2
# q2-mid/right: pairing needs only j<N: -B/sqrt(j)+B/sqrt(N)<0 and -1/(2j)+1/(2N)<0 (j=k or k-1): trivially true; ratio factor decreasing: d/dr (4r+c0)/(r-c0) = -5c0/(r-c0)^2 <0
print("  q2-mid/right: -B/sqrt(j)+B/sqrt(N)<0 and -1/(2j)+1/(2N)<0 as j<N; (4r+c0)/(r-c0)=4+5c0/(r-c0) decreasing; r_a increasing. ANALYTIC")
# mid (i in [2,12]): exponent rate min at i=2: sqrt(3 pi)/(2 sqrt(k-2)) > B/sqrt(N)+0.05/sqrt(N)
A2 = (3*PI).sqrt()/2
lhs2 = A2**2; rhs2 = (B+arb('0.05'))**2
print(f"  mid: A^2=3pi/4={lhs2.str(8)} > (B+0.05)^2={rhs2.str(8)}; A^2(k+1)>(B+0.05)^2(k-2): linear coeff {(lhs2-rhs2).str(6)}>0, at k=99 margin {(lhs2*100-rhs2*97).str(8)}>0 => all k>=99")
assert lhs2 > rhs2 and lhs2*100 - rhs2*97 > 0
# mid2: sqrt(14 pi)/(2 sqrt(k-13)) > (B + 0.05 + 0.10051)/sqrt(N)
A3 = (14*PI).sqrt()/2; C3 = B + arb('0.05') + arb('0.10051')
print(f"  mid2: A^2={ (A3**2).str(8)} > (B+0.15051)^2={(C3**2).str(8)}; A^2(k+1) > C^2(k-13): linear coeff {(A3**2-C3**2).str(6)}>0 => all k>=99")
assert A3**2 > C3**2
# tail: 1/2 > B/sqrt(N) + 1/(2N) + (2k+1)/(2k(k+1)); RHS decreasing, value at k=99:
R99 = B/10 + arb(1)/200 + arb(199)/(2*99*100)
print(f"  tail: RHS(k=99) = {R99.str(10)} < 1/2; RHS decreasing for k>=99 => holds all k>=99")
assert R99 < arb(1)/2
# T2: exponent rate (5/4) sqrt(3 pi)/(2 sqrt k) > (B+0.05)/sqrt(N)
A4 = arb(5)/4*(3*PI).sqrt()/2
print(f"  T2-left: A^2={(A4**2).str(8)} > (B+0.05)^2={((B+arb('0.05'))**2).str(8)} => all k; T2-mid/right: sqrt(3pi)>sqrt(2pi), same pairing as q2. ANALYTIC")
assert A4**2 > (B+arb('0.05'))**2
# T3J: log-deriv <= -ln2 + B/sqrt(N) + 1/(2N) at k=99: -0.693+0.1303 <0, positives decreasing
v = -arb(2).log() + B/10 + arb(1)/200
print(f"  T3J: -ln2 + B/10 + 1/200 = {v.str(10)} < 0 => all k>=99")
assert v < 0
print("Part B: all analytic closures certified. R8 PROVED for all k>=99 (grid redundant).")
