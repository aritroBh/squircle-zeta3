# Exact b_n''(2) for n = 4t+1 and n = 4t+3 from the second-order ODE (see result.md, Sec. 1).
import sympy as sp, json, pickle
from symsum import *
from mpmath import mp, mpf, harmonic, zeta as mz, log as mlog, euler, pi as mpi
B1 = (LN2 - 1)/2                       # b_1'(2)
def b1c(m): return 2/((m+1)**2*(m-1))  # m = 1 mod 4, m >= 5
def b3c(m): return -2/((m-1)**2*(m+1)) # m = 3 mod 4, m >= 3
res = {}
for e in (1, 3):
    n = 4*t + e
    bS = b1c if e == 1 else b3c      # same class as n
    bO = b3c if e == 1 else b1c      # opposite class
    k = 2*i + 1
    T1a = sum_inf(bO(n + 2*k)/k)
    if e == 1:
        T1b = sum_fin(b3c(4*(t-i)-1)/k, 0, t-1) - B1/(2*t+1) - sum_inf(b1c(4*(i-t)+1)/k, t+1)
        O = Hsym(1, 2) - Hsym(1, 1)/2
        T3 = -8*O/(n-1); T4 = sp.Rational(4)/(n+1)
    else:
        T1b = sum_fin(b1c(4*(t-i)+1)/k, 0, t-1) + B1/(2*t+1) - sum_inf(b3c(4*(i-t)-1)/k, t+1)
        O = H_shift(1, 2, 2) - H_shift(1, 1, 1)/2
        T3 = 8*O/(n+1); T4 = sp.Rational(-4)/(n-1)
    T1 = 2*(T1a + T1b)
    m_even = n + 4*i + 4; m_odd = n + 4*i + 2
    T2 = 2*n*bS(n) + 4*(sum_inf(m_even*bS(m_even)) - sum_inf(m_odd*bO(m_odd)))
    b2 = sp.expand((T1 + T2 + T3 + T4)/(1 - n**2))
    syms = sorted(b2.free_symbols - {t}, key=str)
    print('class', e, 'symbols:', syms)
    coll = sp.collect(b2, syms, evaluate=False)
    out = {}
    for key, val in coll.items():
        out[str(key)] = sp.factor(sp.simplify(val))
        print('  [', key, ']:', out[str(key)])
    res[e] = b2
pickle.dump(res, open('b2_law.pkl', 'wb'))

# numeric check against AV8 quadrature data (50-55 digits)
mp.dps = 60
d = json.load(open('../AV8_squircle_sine_fourier/scratchpad/bnder60.json'))
def num(expr, tv):
    sub = {t: tv, LN2: mlog(2), GAM: euler}
    for s_ in range(2, 8): sub[ZETA[s_]] = mz(s_)
    for (s_, a), sym in list(symsum._H.items()):
        sub[sym] = sum(mpf(1)/mpf(j)**s_ for j in range(1, a*tv+1))
    return sp.N(expr.subs(sub), 60)
import symsum
for nn in range(3, 18, 2):
    e = nn % 4; tv = (nn - e)//4
    v = num(res[e], tv)
    print(nn, sp.N(v, 25), 'diff vs AV8:', sp.N(v - sp.Float(d['%d,2' % nn], 60), 3))
