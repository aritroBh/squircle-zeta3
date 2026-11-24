
import json
from mpmath import mp, mpf, zeta, fabs, log, pi, matrix, lu_solve
mp.dps = 135
exec(open('av2core.py').read())

# ---- 1. exact phi_k via sympy series ----
import sympy as sp
w = sp.symbols('w')
gam = sp.Function('gamma')
N = 14
# logGamma(1+w) = -EulerGamma w + sum_{k>=2} (-1)^k zeta(k) w^k/k
lg1 = -sp.EulerGamma*w + sum((-1)**k*sp.zeta(k)*w**k/k for k in range(2, N+2))
lg1 = sp.series(lg1, w, 0, N+1).removeO().expand()
# w^3 G(w)^2/G(2w+1) = w * exp(2*lg1(w) - lg1(2w) - 2*log(w)+2log(w)-... )
# 2 log G(w) - log G(2w+1) = -2 log w + 2 lg1(w) - lg1(2w)
E = sp.series(sp.exp(2*lg1 - lg1.subs(w, 2*w)), w, 0, N+1).removeO().expand()
# psi(2w+1)-psi(w) = 1/w + sum_{k>=2} (-1)^k zeta(k) (2^{k-1}-1) w^{k-1}
H = 1 + sum((-1)**k*sp.zeta(k)*(2**(k-1)-1)*w**k for k in range(2, N+1))
H = sp.series(H, w, 0, N+1).removeO().expand()
phi_ser = sp.series(E*H, w, 0, N+1).removeO().expand()
phi_coeffs = {}
print("=== phi(w) Taylor coefficients (exact) ===")
for k in range(0, N+1):
    c = sp.simplify(phi_ser.coeff(w, k))
    phi_coeffs[k] = c
    print("phi_%d = %s" % (k, c))
mp.dps = 135
# numeric values
z = {k: zeta(k) for k in range(2, 20)}
zneg = {0: mpf(-1)/2, -1: mpf(-1)/12, -2: mpf(0), -3: mpf(1)/120, -4: mpf(0),
        -5: mpf(-1)/252, -6: mpf(0), -7: mpf(1)/240, -8: mpf(0), -9: mpf(-1)/132}
phi_num = {k: mpf(sp.N(phi_coeffs[k], 140)) for k in phi_coeffs}
pred = {  # predicted B2 expansion coefficients c_k = phi_k * zeta(3-k)
    2: None,  # C2 special
    3: phi_num[3]*zneg[0],
    4: phi_num[4]*zneg[-1],
    5: phi_num[5]*zneg[-2],
    6: phi_num[6]*zneg[-3],
    7: phi_num[7]*zneg[-4],
    8: phi_num[8]*zneg[-5],
    9: phi_num[9]*zneg[-6],
    10: phi_num[10]*zneg[-7],
    11: phi_num[11]*zneg[-8],
    12: phi_num[12]*zneg[-9],
}
print("predicted c3 = phi3*zeta(0) =", mp.nstr(pred[3], 30), " (= zeta(3)/2)")
print("predicted c4 = phi4*zeta(-1) =", mp.nstr(pred[4], 30), " (= -3 zeta(4)/16:", mp.nstr(-mpf(3)/16*zeta(4),30),")")
print("predicted c6 =", mp.nstr(pred[6], 30))
print("predicted c8 =", mp.nstr(pred[8], 30))
print("predicted c10 =", mp.nstr(pred[10], 30))
print("predicted c12 =", mp.nstr(pred[12], 30))

# ---- 2. C2 direct at 135 dps ----
C2 = C2_direct(W=mpf(300))
print("C2 direct =", mp.nstr(C2, 60))
print("C2 + zeta(2)/2 =", mp.nstr(C2 + zeta(2)/2, 8))

# ---- 3. B2 at p = 2^k, fit expansion ----
print("=== computing B2 at p=2^k ===")
data = {}
for k in range(8, 25):
    p = 2**k
    q = B2_quad2(p)
    data[k] = q
    print("k=%2d B2-zeta3 = %s" % (k, mp.nstr(q - zeta(3), 12)))
json.dump({str(k): mp.nstr(v, 100) for k, v in data.items()}, open('b2_largep.json','w'))

# fit: R = B2 - zeta3 = c2 e^2 + c3 e^3 + c4 e^4 + c5 e^5 + c6 e^6 + c7 e^7 + c8 e^8 + c9 e^9 + c10 e^10
ks = list(range(14, 25))
unk = [2,3,4,5,6,7,8,9,10]
A = []; b = []
for k in ks:
    e = mpf(1)/2**k
    A.append([e**j for j in unk])
    b.append(data[k] - zeta(3))
Am = matrix(A); bm = matrix(b)
# weighted: multiply row k by 2^{10.5 k} so eps^10 terms ~ O(1)? instead just solve LS
sol = lu_solve(Am.T*Am, Am.T*bm)
resid = Am*sol - bm
print("fit residual max =", mp.nstr(max(fabs(r) for r in resid), 5))
names = {j: "c%d" % j for j in unk}
for j, s in zip(unk, sol):
    print(names[j], "=", mp.nstr(s, 40))
print("--- comparisons ---")
print("c2 vs -zeta(2)/2   :", mp.nstr(sol[0] + zeta(2)/2, 6))
print("c3 vs zeta(3)/2    :", mp.nstr(sol[1] - zeta(3)/2, 6))
print("c4 vs -3zeta(4)/16 :", mp.nstr(sol[2] + mpf(3)/16*zeta(4), 6))
print("c5 (expect 0)      :", mp.nstr(sol[3], 6))
print("c6 vs pred         :", mp.nstr(sol[4] - pred[6], 6), " pred=", mp.nstr(pred[6], 25))
print("c7 (expect 0)      :", mp.nstr(sol[5], 6))
print("c8 vs pred         :", mp.nstr(sol[6] - pred[8], 6), " pred=", mp.nstr(pred[8], 25))
print("c9 (expect 0)      :", mp.nstr(sol[7], 6))
print("c10 vs pred        :", mp.nstr(sol[8] - pred[10], 6), " pred=", mp.nstr(pred[10], 25))
