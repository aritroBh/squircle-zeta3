# ======== STEP 1 verification: exact symbolic f_n, a_n vs 140-digit numeric ========
import sympy as sp, json
from mpmath import mp, mpf, zeta as mzeta
mp.dps = 100
N = 12
eta, s, mu, d = sp.symbols('eta s mu d'); gam = sp.symbols('gam')
Z = sp.symbols(' '.join(f'Z{n}' for n in range(2, 15)))
def zet(n): return Z[n-2]
def T(poly, var, N):
    p = sp.expand(poly); return sum(p.coeff(var, n)*var**n for n in range(0, N+1))
def series_exp(Lc, N):
    Hc = [sp.Integer(1)]
    for m in range(0, N): Hc.append(sp.expand(sum((j+1)*Lc[j+1]*Hc[m-j] for j in range(0, m+1))/(m+1)))
    return Hc
Lc = [sp.Integer(0)]*(N+1); Lc[1] = -2*gam*mu
for n in range(2, N+1):
    for j in range(0, n+1, 2): Lc[n] += 2*(-1)**n*zet(n)/n*sp.binomial(n, j)*mu**(n-j)*s**j
Htot = sum(series_exp([sp.expand(c) for c in Lc], N))
cc = [sp.Integer(1)]
for k in range(1, N//2+1): cc.append(sp.expand((1/sp.factorial(2*k) - sum(cc[i]*cc[k-i] for i in range(1, k)))/2))
logKc = [sp.Integer(0)]*(N+1); logKc[1] = gam
for n in range(2, N+1): logKc[n] = -(-1)**n*zet(n)/n
tmp = series_exp(logKc, N)
invG2 = [sp.expand(sum(tmp[j]*(-1)**(m-j) for j in range(0, m+1))) for m in range(N+1)]
main = [sp.Integer(0)]*(N+1)
for k in range(0, N//2+1):
    H2k = sp.expand((sp.factorial(2*k)*sp.expand(Htot).coeff(s, 2*k)).subs(mu, eta/2))
    for m in range(0, N+1-2*k): main[2*k+m] += cc[k]*sp.expand(H2k).coeff(eta, m)
Fc = [sp.expand(sum(invG2[j]*main[m-j] for j in range(0, m+1))) for m in range(N+1)]
eta_ser = sum((-d)**j for j in range(1, 13)); eta_pow = [sp.Integer(1)]
for m in range(1, 12): eta_pow.append(T(eta_pow[-1]*eta_ser, d, 11))
Ac_eta = [sp.expand(4*(Fc[m] + (Fc[m-1] if m >= 1 else 0))) for m in range(0, 11)]
Ac = [sp.expand(sum(Ac_eta[m]*eta_pow[m] for m in range(0, 11))).coeff(d, n) for n in range(0, 11)]

zval = {zet(n): sp.Float(mp.nstr(mzeta(n), 90), 90) for n in range(2, 13)}
num = json.load(open('step2_num.json'))
print("EXACT(symbolic) - NUMERIC(140dps engine), should be ~1e-88:")
worst = 0
for n in range(0, 11):
    fv = mpf(str(Fc[n].subs(zval).evalf(90))); fn = mpf(num['f'][n])
    av = mpf(str(Ac[n].subs(zval).evalf(90))); an = mpf(num['a'][n])
    df, da = abs(fv-fn)/max(1,abs(fn)), abs(av-an)/max(1,abs(an))
    worst = max(worst, df, da)
    print(f"n={n:2d}  rel diff f: {mp.nstr(df,3)}   rel diff a: {mp.nstr(da,3)}")
print("worst relative residual:", mp.nstr(worst, 3))
