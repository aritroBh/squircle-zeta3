import sympy as sp
# Asymptotic coefficients h_r of H_m = (-1)^m sum_{r>=2} h_r (m pi)^{-r}
v, m = sp.symbols('v m', positive=True)
pmax = 26
def CA(j): return sp.integrate(v**j*sp.cos(v), (v, 0, sp.pi))
def CB(j): return sp.integrate(v**j*sp.sin(v), (v, 0, sp.pi))
CAv = [sp.simplify(CA(j)) for j in range(pmax)]
CBv = [sp.simplify(CB(j)) for j in range(pmax)]
# H_m = -(-1)^m sum_k (-1)^k [(2k)! A_{2k+2} + (2k+1)! B_{2k+3}]
# A_p = sum_j (-1)^j binom(p+j-1,j) (mπ)^{-p-j} CA_j ; B_p similarly with CB_j, sign (-1)^j as well
KMAX = 12
h = {}
for k in range(KMAX):
    for (P, moments, fk) in [(2*k+2, CAv, sp.factorial(2*k)), (2*k+3, CBv, sp.factorial(2*k+1))]:
        for j in range(0, pmax - P):
            r = P + j
            if r > pmax: continue
            coeff = (-1)**j * sp.binomial(P+j-1, j) * moments[j]
            h[r] = h.get(r, 0) + (-1)*(-1)**k * fk * coeff
h = {r: sp.simplify(h[r]) for r in sorted(h)}
for r in sorted(h)[:14]:
    print(f"h_{r} = {h[r]}")
import json
json.dump({str(r): str(h[r]) for r in h}, open("h_coeffs.json","w"))
print("saved h_coeffs.json")
