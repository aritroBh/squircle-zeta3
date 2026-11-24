import sympy as sp, json
pmax = 30
CA = [0]*pmax; CB = [0]*pmax
CA[0] = sp.Integer(0); CB[0] = sp.Integer(2)
for j in range(1, pmax):
    CA[j] = sp.simplify(-j*CB[j-1])
    CB[j] = sp.simplify(sp.pi**j + j*CA[j-1])
KMAX = 12
h = {}
for k in range(KMAX):
    for (P, moments, fk) in [(2*k+2, CA, sp.factorial(2*k)), (2*k+3, CB, sp.factorial(2*k+1))]:
        for j in range(0, pmax - P):
            r = P + j
            if r > pmax: continue
            coeff = (-1)**j * sp.binomial(P+j-1, j) * moments[j]
            h[r] = h.get(r, 0) + (-1)*(-1)**k * fk * coeff
h = {r: sp.simplify(h[r]) for r in sorted(h)}
for r in sorted(h):
    print(f"h_{r} = {sp.nsimplify(h[r])}   ~{sp.N(h[r],8)}")
json.dump({str(r): str(h[r]) for r in h}, open("h_coeffs.json","w"))
print("saved")
