code = r'''# step_04.py: the zeta(5) edge law.  Z_e^(5) = 5!*[z5] int E5_local = 120*24*int G,  G = g41+g32+g23+g14
# (because [z5] int t^n X^a Y^b dt = 24/(n+1) = 24 int t^n for a+b=5, a,b>=1, and all other local
# monomials have weight <= 4).  Z_e^(2,3) = 5!*[z2z3] = 120*(-12)*int(g32+g23).
import sympy as sp, pickle

t, al, be, L, h, s0, u, v = sp.symbols('t alpha beta L h s0 u v', real=True)
R = pickle.load(open('steps/step_01_reduced.pkl', 'rb'))
sy = {'t': t, 'alpha': al, 'beta': be, 'L': L, 'h': h, 's0': s0}
Ig = {k: sp.sympify(v, locals=sy) for k, v in R['Ig'].items()}
G41, G32, G23, G14 = Ig['(4, 1)'], Ig['(3, 2)'], Ig['(2, 3)'], Ig['(1, 4)']

Z5 = sp.expand(2880*(G41 + G32 + G23 + G14))
Z23 = sp.expand(-1440*(G32 + G23))
print("int G =", sp.factor(G41 + G32 + G23 + G14))
print("int(g32+g23) =", sp.factor(G32 + G23))

# common structure: everything over L^4? find denominator
den5 = sp.denom(sp.together(Z5)); den23 = sp.denom(sp.together(Z23))
print("denom Z5:", den5, " denom Z23:", den23)

# Re-express with u = -s0/L, v = 1 + s0/L (i.e. s0 = -uL, L+s0 = vL) and h^2/L pulled out.
# substitute s0 -> -u*L, then L + s0 -> L*v automatically becomes L(1-u); handle via v = 1-u.
Z5u = sp.expand(Z5.subs(s0, -u*L))
Z23u = sp.expand(Z23.subs(s0, -u*L))
# replace (1-u) by v to symmetrize
Z5uv = sp.expand(Z5u.subs(u, 1 - v)).subs(v, 1 - u)  # identity; instead do manual below
# Manual symmetric form: write everything as polynomial in u and (1-u), then set v = 1-u.
Z5uv = sp.expand(Z5u)
Z23uv = sp.expand(Z23u)
# factor out alpha*beta*h^2/L^? : determine overall factor by dividing
f5 = sp.simplify(Z5uv/(al*be*h**2/L**4))
f23 = sp.simplify(Z23uv/(al*be*h**2/L**4))
print("Z5/(alpha*beta*h^2/L^4) =", sp.factor(f5))
print("Z23/(alpha*beta*h^2/L^4) =", sp.factor(f23))
f5L = sp.simplify(Z5uv/(al*be*h**2/L))
f23L = sp.simplify(Z23uv/(al*be*h**2/L))
print("Z5/(alpha*beta*h^2/L) =", sp.factor(f5L))
print("Z23/(alpha*beta*h^2/L) =", sp.factor(f23L))
with open('steps/step_04_out.txt', 'w') as f:
    f.write(f"int G = {sp.factor(G41+G32+G23+G14)}\nint(g32+g23) = {sp.factor(G32+G23)}\n")
    f.write(f"Z5 raw = {Z5}\nZ23 raw = {Z23}\n")
    f.write(f"Z5/(ab h^2/L) = {sp.factor(f5L)}\nZ23/(ab h^2/L) = {sp.factor(f23L)}\n")
print("saved steps/step_04_out.txt")
'''
open('steps/step_04.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_04.py'],capture_output=True,text=True,timeout=3000)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
