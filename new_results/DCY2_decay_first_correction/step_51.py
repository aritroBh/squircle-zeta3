import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
audit = r'''"""Audit every analytic bound against the true value (mpmath, high precision)."""
from mpmath import mp, mpf, pi, sqrt, exp, gamma, log, beta, gammainc
j = 1j
omega = exp(j*pi/4); wb = exp(-j*pi/4)
bR = pi/(2*sqrt(2)); b = pi*wb/2
def Hc(s): return beta(s/2, mpf(3)/2)/2
def df2(k):
    r = mpf(1)
    while k > 1: r *= k; k -= 2
    return r
kap0 = 2.48

def true_quant(N):
    N = mpf(N); q = {}
    q['y'] = y = sqrt(pi/(2*N)); q['rho'] = rho = 1/y; q['s0'] = s0 = rho*omega
    q['theta'] = th = sqrt(pi*N)
    q['kap'] = kap = sqrt(kap0**2 + 4*log(N/400))
    def ell(m):
        if m == 1: return (N-1)*y - N*y/(1+y*wb) + j*pi*omega/2
        return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
    q['l'] = l = {m: ell(m) for m in range(1,7)}
    l1,l2 = l[1],l[2]
    q['c2'] = -l2.real
    q['c2min'] = c2m = 0.9*bR*y
    q['sig'] = sig = 1/sqrt(2*c2m)
    q['Delta'] = Delta = kap*sig
    def ellt(t): return (N-1)*log(1+y*t) - N*log(1+y*t/(1+y*wb)) + j*pi*omega*t/2
    q['ellt'] = ellt
    def R3(t): return ellt(t) - l1*t - l2*t**2
    q['R3'] = R3
    q['F'] = lambda t: Hc(s0+omega*t)/Hc(s0)
    # true |cdot| moments of the Gaussian
    def Mtrue(jj):
        return mp.quad(lambda t: abs(t)**jj*exp(l1.real*t + l2.real*t**2), [-mp.inf, 0, mp.inf])
    q['Mtrue'] = Mtrue
    def TPtrue(p):
        return (mp.quad(lambda t: t**p*exp(l1.real*t+l2.real*t**2), [Delta, mp.inf])
                + mp.quad(lambda t: (-t)**p*exp(-l1.real*t+l2.real*t**2), [Delta, mp.inf]))
    q['TPtrue'] = TPtrue
    # true P polynomial and E
    def Hder(m):
        return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
    H0 = Hc(s0)
    h = [1] + [omega**m*Hder(m)/H0/mp.factorial(m) for m in range(1,4)]
    l3,l4,l5 = l[3],l[4],l[5]
    P = {}
    def add(p,c): P[p] = P.get(p,0)+c
    add(0,1); add(1,h[1]); add(2,h[2]); add(3,h[3])
    add(3,l3); add(4,l4); add(5,l5)
    add(6,l3**2/2); add(7,l3*l4); add(9,l3**3/6)
    add(4,h[1]*l3); add(5,h[1]*l4); add(5,h[2]*l3); add(7,h[1]*l3**2/2)
    q['P'] = P
    q['E'] = lambda t: q['F'](t)*exp(R3(t)) - sum(cp*t**p for p,cp in P.items())
    # true mu_p^inf and Model
    t0 = -l1/(2*l2)
    def mu(p):
        base = exp(-l2*t0**2)*sqrt(-pi/l2)
        return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
    q['Model'] = (2/pi)*exp((N-1)*log(s0)-N*log(1+s0)+j*pi*s0/2)*H0*omega*sum(cp*mu(p) for p,cp in P.items())
    q['It'] = (2/pi)*mp.quad(lambda x: Hc(omega*x)*exp(j*pi*omega*x/2)*(omega*x/(1+omega*x))**N/(omega*x)*omega,
                             [0, rho/8, rho/2, rho, 2*rho, 4*rho, 10*rho, mp.inf])
    # true tails
    def Lam(t):
        u = 1+y*t
        Q = u**2+sqrt(2)*y*u+y**2
        return (N-1)*log(u) + (N/2)*log(1+sqrt(2)*y+y**2) - (N/2)*log(Q) - bR*t
    q['Lam'] = Lam
    q['TR_true'] = mp.quad(lambda t: abs(Hc(s0+omega*t))*exp(Lam(t)), [Delta, mp.inf])*exp(((N-1)*log(s0)-N*log(1+s0)+j*pi*s0/2).real)
    return q

def audit(N):
    import dcy2_budget as B
    B.P.kap0 = kap0
    qq = B.quant(N); tq = true_quant(N)
    out = []
    # 1. E bound: sample |E(t)| vs B_j |t|^j bound: recompute the coefficient list
    # (reuse: E_mom vs true int |E| e^{Re} dt)
    Eint_true = mp.quad(lambda t: abs(tq['E'](t))*exp(tq['l'][1].real*t+tq['l'][2].real*t**2),
                        [-tq['Delta'], 0, tq['Delta']])
    out.append(("E_mom >= int|E| e^{Re}", qq['E_mom'], Eint_true))
    # 2. |cdot| moments
    for jj in [0,2,4,6,8,10,12]:
        # recompute Mtilde from budget internals: quick re-implementation
        l1,l2 = tq['l'][1],tq['l'][2]
        c2m = qq['c2min']; ReL1 = abs(l1.real); tR = ReL1/(2*c2m)
        efac = exp(ReL1**2/(4*c2m))
        Mt = efac*sum(mp.binomial(jj,i)*tR**(jj-i)*c2m**(-(i+1)/2)*gamma((i+1)/2) for i in range(jj+1))
        out.append((f"Mtilde({jj}) >= true", Mt, tq['Mtrue'](jj)))
    # 3. TP_p
    from mpmath import gammainc as gi
    for p in [0,1,2,3,9]:
        c2m = qq['c2min']; Delta = qq['Delta']
        TP = mpf(1)/2*(c2m/2)**(-(p+1)/2)*gi((p+1)/2, c2m*Delta**2/2, regularized=False)
        out.append((f"TP({p}) >= true", TP, tq['TPtrue'](p)))
    # 4. right tail
    out.append(("T_R >= true right tail", qq['Tails'], tq['TR_true']))
    # 5. total
    out.append(("E_total >= |It - Model|", qq['E_total'], abs(tq['It']-tq['Model'])))
    # 6. budget
    out.append(("budget >= E_total", qq['budget'], qq['E_total']))
    print(f"=== N = {N} ===")
    for name, bd, tv in out:
        ok = bd >= tv
        print(f"  [{'OK ' if ok else 'FAIL'}] {name}: bound={float(bd):.3e} truth={float(tv):.3e} ratio={float(tv/bd):.3f}")
    return out

if __name__ == "__main__":
    mp.dps = 50
    for N in [400, 2000, 10000]:
        audit(N)
'''
open("dcy2_audit.py","w").write(audit)
from mpmath import mp
mp.dps = 50
import dcy2_audit
for N in [400, 2000, 10000]:
    dcy2_audit.audit(N)
