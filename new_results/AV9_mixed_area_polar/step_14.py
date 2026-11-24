
# step_14.py : square-end self-dual to order 10 (extend table; test zeta(9) first appearance)
import sympy as sp, pickle, time
import av9core as C
t0=time.time()
N = 10
h, e = C.selfdual_coeffs(N)
tab = C.beta_int_table(N)
out = {}
for n in range(9, N+1):
    out[n] = sp.expand(C.integrate_poly(4*e[n], tab))
    print(f"c_{n} done t={time.time()-t0:.0f}s", flush=True)
    print(f"c_{n} =", out[n])
    print("  numeric:", sp.N(out[n], 30))
    print("  zeta(9) coeff:", out[n].coeff(sp.zeta(9),1), " zeta(7):", out[n].coeff(sp.zeta(7),1), flush=True)
pickle.dump({n: str(out[n]) for n in out}, open('scratchpad/c_exact_9_10.pkl','wb'))
print("TOTAL", time.time()-t0)
