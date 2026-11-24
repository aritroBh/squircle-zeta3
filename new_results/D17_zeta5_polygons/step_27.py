import os, subprocess
d17 = "research/08_theory/D17_zeta5_polygons"
venv_py = "./.venv/bin/python"

script = r'''
import sympy as sp, pickle
phi = sp.Symbol('phi')
h, Sh, Ch = sp.symbols('h Sh Ch', positive=True)
data = pickle.load(open('monolist5.pkl','rb'))
terms = [(b,p,q,sp.sympify(r)) for b,p,q,r in data]
print("The 18 terms (cos^b phi)(Sh^p)(Ch^q)(coeff):")
for b,p,q,r in terms:
    print("  b=%d p=%d q=%d coeff=%s"%(b,p,q,sp.factor(r)))
def INT(b):
    if b==-2: return 2*sp.tan(h)
    if b==-4: return 2*(sp.tan(h)+sp.tan(h)**3/3)
    if b==-6: return 2*(sp.tan(h)+sp.Rational(2,3)*sp.tan(h)**3+sp.tan(h)**5/5)
I = 0
for b,p,q,r in terms:
    I += r*Sh**p*Ch**q*INT(b)
I = sp.expand(I).subs(sp.tan(h), Sh/Ch)
I = sp.expand(I)
# reduce Ch^2 -> 1-Sh^2
Ir = sp.expand(I).subs(Ch**2, 1-Sh**2)
Ir = sp.expand(Ir)
for pw in [3,4,5,6,7,8]:
    Ir = sp.expand(Ir.subs(Ch**pw, Ch**(pw-2)*(1-Sh**2)))
Ir = sp.simplify(Ir)
print("\nJ5(h) = int_{-h}^{h} (G41+G32+G23+G14) dphi =")
print(sp.factor(Ir))
with open('J5_of_h.txt','w') as f: f.write(sp.srepr(sp.factor(Ir)))
# C_m = 5760 m J5, m = pi/(2h)
m,pi = sp.symbols('m pi', positive=True)
Cm = sp.expand(5760*(pi/(2*h))*Ir)
Cm = sp.simplify(Cm)
print("\nC_m^(5) = 5760 m J5  (m=pi/2h) =")
print(sp.factor(Cm))
# Check m=2: h=pi/4, Sh=Ch=sqrt2/2
import mpmath as mp
mp.mp.dps=50
fJ = sp.lambdify((Sh,), Ir.subs(Ch, sp.sqrt(1-Sh**2)), 'mpmath')
J2 = fJ(mp.sin(mp.pi/4))
print("\nJ5(h=pi/4) =", mp.nstr(J2,30), " target -sqrt2/128 =", mp.nstr(-mp.sqrt(2)/128,30))
print("resid:", mp.nstr(abs(J2+mp.sqrt(2)/128),3))
print("C_2 = 11520*J5 =", mp.nstr(11520*J2,20), " target -90sqrt2 =", mp.nstr(-90*mp.sqrt(2),20))
'''.replace("D17DIR", d17)
with open(os.path.join(d17,"integ5.py"),"w") as f: f.write(script)
r = subprocess.run([venv_py, os.path.join(d17,"integ5.py")], capture_output=True, text=True, cwd=d17, timeout=600)
print(r.stdout); print("ERR:", r.stderr[-3000:])
