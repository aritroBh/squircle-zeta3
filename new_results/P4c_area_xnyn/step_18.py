# step_16.py : R3 verification (both d3 candidates) + C(n) O(e) coefficient extraction
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta
import pickle
mp.dps = 70
L2 = log(2); z2 = pi**2/6; z3 = zeta(3)
b3 = mpf('0.5') - L2/2
c3 = L2**2/4 - L2/2 + 1 + z2
d3_prev = -L2**3/24 + L2**2/8 - L2/2 - L2*z2/2 + z2/2 - z3/2 + 1 + 3*sqrt(2)/2   # +3sqrt2/2
d3_mine = d3_prev - sqrt(2)/4                                                   # +5sqrt2/4
print("d3_prev =", nstr(d3_prev, 50)); print("d3_mine =", nstr(d3_mine, 50))
def R3(n, d3):
    ln = log(mpf(n))
    return ln**3/3 + b3*ln**2 + c3*ln + d3
d = pickle.load(open("E_values_new.pkl","rb"))
print()
print("R3_num = n*resid ; diffs with two candidates:")
for k,v in d.items():
    n = mpf(k); E = mpf(v[0]); resid = mpf(v[1])
    R3num = n*resid
    dp = R3num - R3(n, d3_prev); dm = R3num - R3(n, d3_mine)
    print("n=%s: R3num=%s" % (k, nstr(R3num, 40)))
    print("   diff(prev d3)=%s   diff(mine)=%s   (log^4 n)/n = %s" % (nstr(dp,16), nstr(dm,16), nstr(log(n)**4/n,8)))

# C(n) from the runs: C = n^2 (Phi/(1+1/2n) - E) ; compare with T1 + e(T1-T2) [+ delta]
def mtil(k):  # m̃_k = Σ_{i<=k} e_i, from E_g series: closed forms
    pass
# e1=-L2, e2=L2^2/2+z2/2, e3=-z3/4-L2*z2/2-L2^3/6, e4=z4/4+L2*z3/4+L2^2*z2/4+z2^2/8+L2^4/24
z4 = pi**4/90
e1 = -L2; e2 = L2**2/2+z2/2; e3 = -z3/4-L2*z2/2-L2**3/6; e4 = z4/4+L2*z3/4+L2**2*z2/4+z2**2/8+L2**4/24
eks = [mpf(1), e1, e2, e3, e4]
mt = [sum(eks[:k+1]) for k in range(5)]
def T1(N):
    L = log(mpf(N))
    return L**2/2 + mt[1]*L + mt[2]
def T2(N):
    L = log(mpf(N))
    return L**3/3 + mt[1]*L**2 + 2*mt[2]*L + 2*mt[3]
print()
print("C(n) vs T1 + e(T1-T2) + delta*(1/(2sqrt2)):")
for k,v in d.items():
    n = mpf(k); C = mpf(v[2]); e = 1/n; N = n*sqrt(2)
    base = T1(N) + e*(T1(N)-T2(N))
    print("n=%s: C=%s" % (k, nstr(C, 40)))
    print("   C - base = %s   [+e/(2sqrt2)] residual = %s" % (nstr(C-base, 18), nstr(C-base+e/(2*sqrt(2)), 18)))
