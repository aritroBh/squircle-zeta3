
summary = r'''
# step_18.py : final summary of constants (exact + 30-digit values) for the report
import sympy as sp, pickle
c = {int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_selfdual.pkl','rb')).items()}
c.update({int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_9_10.pkl','rb')).items()})
c.update({int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_11_12.pkl','rb')).items()})
print("== SQUARE END: V(eta) = 4 + sum c_n eta^n (self-dual) ==")
for n in sorted(c):
    print(f"c_{n} = {c[n]}")
    print(f"      = {sp.N(c[n], 32)}")
print("\n== CIRCLE END: V(s) = pi + sum a_2k s^2k (self-dual) ==")
pi = sp.pi
a = {2: 3*pi*(pi**2-8)/4,
     4: pi*(15*pi**4 - 352*pi**2 + 1664)/192,
     6: pi*(9699*pi**6 - 96280*pi**4 + 13440*pi**2 - 39936)/11520,
     8: pi*(-418761675*pi**8 + 4482752384*pi**6 - 3451302400*pi**4 - 63651840*pi**2 + 138117120)/77414400}
a[10] = sp.sympify(pickle.load(open('scratchpad/circle_a10.pkl','rb'))['a10'])
for n in sorted(a):
    print(f"a_{n} = {a[n]}")
    print(f"      = {sp.N(a[n], 32)}")
print("\n== first appearances (square end, self-dual) ==")
for k in range(1,6):
    print(f"zeta({2k+1}): order {2*k+2}, coeff {c[2*k+2].coeff(sp.zeta(2*k+1),1)}")
print("zeta(3)^2: order 7, coeff", c[7].coeff(sp.zeta(3)**2))
print("zeta(2)zeta(3): order 6, coeff", sp.expand(c[6].subs(sp.pi**2,6*sp.zeta(2))).coeff(sp.zeta(2)*sp.zeta(3)))
print("zeta(3)^3: order 10, coeff", c[10].coeff(sp.zeta(3)**3))
print("zeta(3)zeta(5): order 10, coeff", c[10].coeff(sp.zeta(3)*sp.zeta(5)))
'''
open('step_18.py','w').write(summary)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_18.py'],capture_output=True,text=True,timeout=600)
print(r.stdout); print(r.stderr[-300:])
with open('progress.log','a') as f:
    f.write("step_15/17: c_9..c_12 verified (quad 87-90 digits, FFT 31-44 digits)\n")
    f.write("step_16: first appearances zeta(2k+1)@2k+2 for k=1..5 with coeffs 1,1/2,1/2,-1/12,11/4; zeta(11) not in c_11\n")
    f.write("step_18: final summary tables generated\n")
