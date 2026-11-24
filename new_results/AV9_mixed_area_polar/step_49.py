
step16 = r'''
# step_16.py : square-end to order 12; first-appearance coefficients of zeta(2k+1)
import sympy as sp, pickle, time
import av9core as C
t0=time.time()
N = 12
h, e = C.selfdual_coeffs(N)
tab = C.beta_int_table(N)
out = {}
for n in range(1, N+1):
    out[n] = sp.expand(C.integrate_poly(4*e[n], tab))
    print(f"c_{n} done t={time.time()-t0:.0f}s", flush=True)
print("\nfirst-appearance coefficients [zeta(2k+1)] c_{2k+2}:")
for k in range(1, 6):
    n = 2*k+2
    cf = out[n].coeff(sp.zeta(2*k+1),1)
    print(f"  [zeta({2k+1})] c_{n} = {cf}")
print("\n[zeta(11)] c_11 =", out[11].coeff(sp.zeta(11),1), " (should be 0)")
print("[zeta(3)^3] c_10 =", out[10].coeff(sp.zeta(3)**3))
print("c_11 numeric:", sp.N(out[11],25))
print("c_12 numeric:", sp.N(out[12],25))
pickle.dump({n: str(out[n]) for n in [11,12]}, open('scratchpad/c_exact_11_12.pkl','wb'))
print("TOTAL", time.time()-t0)
'''
open('step_16.py','w').write(step16)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_16.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout); print(r.stderr[-400:])
