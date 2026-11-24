# claude_check_lead2.py -- independent numeric checks of LEAD2's new identities (M1 session, 2026-10-27).
# Own code (not LEAD2's lib): Master Lemma (i),(iii); Lemma 5.2 (H = B(1+x,y)); Theorem 5.4 (T_d);
# (7) (u m_e part); the end-to-end K_d = (d!/2) e_{d+1} from the u-measures, d = 1..7; and the
# measure dictionary (pole combinations in v -> u-measures) at random points.
import mpmath as mp
mp.mp.dps = 40
E = lambda s: mp.gamma(1 + s/2)*mp.gamma(mp.mpf(1)/2 - s/2)/mp.sqrt(mp.pi)
def e_coef(n):  # [y^n] E(y) by Taylor
    return mp.taylor(E, 0, n)[n]
out = []
def rep(name, r, tol=mp.mpf(10)**-25):
    ok = abs(r) < tol; out.append(ok); print(('PASS' if ok else 'FAIL'), name, mp.nstr(r, 3), flush=True)

q = lambda w: (1 + w)/(1 + w*w)**mp.mpf(1.5)
h = lambda w: ((1 + w)/mp.sqrt(1 + w*w) - 1)/w
g1 = lambda w: (1 - w)/(1 + w)*h(w)
for s in [mp.mpf('0.37'), mp.mpc('-0.2', '0.3')]:
    WW = mp.quad(lambda w: w**(-s)*q(w), [0, 1, mp.inf])
    GG = mp.quad(lambda w: w**(-s)*g1(w), [0, 1, mp.inf])
    rep(f'WW({s}) = E(s)+E(-s)', WW - (E(s) + E(-s)))
    rep(f'GG({s}) = 2pi/sin(pi s) - (E(s)+E(-s))/s', GG - (2*mp.pi/mp.sin(mp.pi*s) - (E(s) + E(-s))/s))
x, y = mp.mpf('0.3'), mp.mpf('0.4')
H = mp.quad(lambda u: ((1 + u)/2)**x*((1 - u)/2)**y/(1 - u), [-1, 0, 1])
rep('H(0.3,0.4) = B(1.3,0.4)', H - mp.beta(1 + x, y))

L1 = lambda u: mp.log((1 + u)/2); L2 = lambda u: mp.log((1 - u)/2)
S = lambda u: mp.sqrt(1 + u*u)
m_e = lambda u: (1 - mp.sqrt(2)/S(u))/(1 - u*u) if abs(1 - u*u) > mp.mpf(10)**-30 else mp.mpf(-1)/4
m_odd = lambda u: (1 + mp.sqrt(2)*u/S(u))/(1 - u*u)
m_diag = lambda u: (1 - u)/2*m_e(u)
m_oddlast = lambda u: m_odd(u) - 1/(2*(1 - u))
M = lambda i, j, m: mp.quad(lambda u: L1(u)**i*L2(u)**j*m(u), [-1, 0, 1])
A = lambda p, qq: M(p, qq, lambda u: 1/(1 - u))
for d in [1, 3, 5]:
    mid = (d - 1)//2
    T = sum((-1)**(j + 1)*mp.binomial(d, j)*A(d - j, j) for j in range(mid + 1)) - mp.mpf(1)/2*(-1)**(mid + 1)*mp.binomial(d, mid)*A(mid + 1, mid)
    rep(f'T_{d}/d! = (1-2^-d) zeta(d+1)', T/mp.factorial(d) - (1 - mp.mpf(2)**-d)*mp.zeta(d + 1))
    um = sum((-1)**(j + 1)*mp.binomial(d, j)*M(d - j, j, lambda u: u*m_e(u)) for j in range(mid + 1))
    rep(f'(u m_e part)_{d}/d! = (1-2^-d)zeta(d+1) - e_(d+1)/2', um/mp.factorial(d) - ((1 - mp.mpf(2)**-d)*mp.zeta(d + 1) - e_coef(d + 1)/2))
for d in range(1, 8):
    K = 0
    for j in range(d//2 + 1):
        if d % 2 == 0:
            m = m_diag if j == d//2 else m_e
        else:
            m = m_oddlast if j == (d - 1)//2 else m_odd
        K += (-1)**(j + 1)*mp.binomial(d, j)*M(d - j, j, m)
    rep(f'K_{d} (eq. (3), diagonal included) = (d!/2) e_(d+1)', K - mp.factorial(d)/2*e_coef(d + 1))
# measure dictionary: (pole combination in v) * dv/du = m(u) du at random u
a, b = 1 + mp.sqrt(2), mp.sqrt(2) - 1
for u in [mp.mpf('-0.71'), mp.mpf('0.13'), mp.mpf('0.88')]:
    v = u + S(u); jac = v/S(u)
    rep(f'm_odd at u={u}', (1/(v + a) - 1/(v - a))*jac - m_odd(u))
    rep(f'm_e at u={u}', (1/(v + a) - 1/(v + b))*jac - m_e(u))
    rep(f'm_diag at u={u}', (1/(v + a) - 1/(2*v))*jac - m_diag(u))
    rep(f'm_oddlast at u={u}', (1/(v + a) - 1/(2*(v - a)) + 1/(2*(v + b)) - 1/(2*v))*jac - m_oddlast(u))
    rep(f't = (1+u)/2 at u={u}', (v*v + 2*v - 1)/(4*v) - (1 + u)/2)
print(f'{sum(out)}/{len(out)} PASS', 'ALL PASS' if all(out) else 'FAILURES')
