# Numerical re-verification (>= 70 digits) of (i) every depth-2 reduction from dshuffle.py, including the
# classical evaluation of the one free weight-4 element (not needed for R, checked anyway), (ii) exact R against
# SQSR's independent 260-digit Euler-Maclaurin value.  Z computed from an integral representation (quad).
import sympy as sp, pickle, time
from mpmath import mp, mpf, quad, gamma, zeta, log, pi, polylog, nstr
from dshuffle import reduce_weight, ZS, l2
mp.dps = 90
def Znum(a, s, b, t):
    # Z(a,s;b,t) = s/Gamma(a) int_0^1 (-log u)^(a-1) Li_b(s t u)/(1 - s u) du   (from 1/m^a = int u^(m-1)(-log u)^(a-1)/Gamma(a))
    f = lambda u: (-log(u))**(a-1)*polylog(b, s*t*u)/(1-s*u)
    return s/gamma(a)*quad(f, [0, mpf(1)/2, 1])
consts = {l2: log(2)}
for s_ in range(2, 8): consts[ZS[s_]] = zeta(s_)
worst = 0
t0 = time.time()
for w in range(2, 6):
    red, free = reduce_weight(w)
    for sym, val in red.items():
        nm = str(sym)[2:]; a, s, b, t = int(nm[0]), (1 if nm[1] == 'p' else -1), int(nm[3]), (1 if nm[4] == 'p' else -1)
        x = Znum(a, s, b, t)
        if val.free_symbols & free:
            # classical: Z(3,-;1,-) = sum_{m>l} (-1)^{m+l}/(m^3 l); from sum (-1)^{n+1}H_n/n^3 =
            # 11/4 z4 - 7/4 z3 log2 + pi^2/12 log^2 2 - log^4 2/12 - 2 Li4(1/2)
            S = mpf(11)/4*zeta(4) - mpf(7)/4*zeta(3)*log(2) + pi**2/12*log(2)**2 - log(2)**4/12 - 2*polylog(4, mpf(1)/2)
            z31 = Znum(3, -1, 1, -1)
            # sum (-1)^{n+1} H_n/n^3 = sum (-1)^{n+1}/n^3 (sum_{l<n} 1/l) + eta(4) = -Z(3,-;1,+) + eta(4)
            chk = -Znum(3, -1, 1, 1) + (1-mpf(2)**-3)*zeta(4) - S
            if str(sym) == 'Z_3m_1m': print('   classical sum (-1)^{n+1}H_n/n^3 formula: residual', nstr(chk, 3))
            y = sp.N(val.subs({**consts}).subs(sym, 0), 95)
            print('  ', sym, '(free; relation checked with it numerically)')
            continue
        y = mpf(str(sp.N(val.subs(consts), 95)))
        d = abs(x-y); worst = max(worst, d)
        print('  w=%d' % w, sym, nstr(x, 30), ' residual', nstr(x-y, 3), flush=True)
# check the weight-4 relations containing the free symbol, with its numerical value
red4, free4 = reduce_weight(4)
fz = list(free4)[0]; fv = Znum(3, -1, 1, -1)
for sym, val in red4.items():
    nm = str(sym)[2:]; a, s, b, t = int(nm[0]), (1 if nm[1] == 'p' else -1), int(nm[3]), (1 if nm[4] == 'p' else -1)
    x = Znum(a, s, b, t); y = mpf(str(sp.N(val.subs(consts).subs(fz, sp.Float(str(fv), 95)), 95)))
    worst = max(worst, abs(x-y))
    print('  w=4 (with free value)', sym, ' residual', nstr(x-y, 3))
print('worst Z residual', nstr(worst, 3), ' time %.0fs' % (time.time()-t0))
ex = pickle.load(open('r3_R_exact.pkl', 'rb'))
mp.dps = 270
Rsq = mpf('-0.0092253434023847181165512739452309923732365777972675358869518389291408934672677754473923510209906214055720303570532840333831215324734087394722327881111739581685694818145258345004043867302817013857584807010239733411379671548439807706582233658184483293188223947629')
mp.dps = 270
yR = mpf(str(sp.N(ex['R'].subs({l2: sp.log(2), **{ZS[s]: sp.zeta(s) for s in range(2, 8)}}), 270)))
print('exact R vs SQSR 260-digit Euler-Maclaurin value: diff', nstr(yR-Rsq, 3))
