import sympy as sp
t=sp.symbols('t', positive=True); r2=sp.sqrt(2)
expr=(t**2+r2*t)**2+(r2*t+1)**2-(t**2+r2*t+1)**2
print("p2mag modulus identity residual (sympy expand):", sp.expand(expr))
# also verify the H2 algebra: [2(x+rt2 k)]^2 - (x^2+2 rt2 k x+4k^2)
x,k=sp.symbols('x k', positive=True)
print("H2 2Re w >= |w| numerator:", sp.expand((2*(x+r2*k))**2-(x**2+2*r2*k*x+4*k**2)))
# AM-GM residual
print("AM-GM: 1+r2 x+x^2-(2+r2)x =", sp.expand(1+r2*x+x**2-(2+r2)*x))
with open('progress.log','a') as f:
    f.write("turn5 R2: grep+diff audit done. Only intended fix lines changed vs buggy baselines (step_18 unchanged).\n"
            "  ecw fixed = 0.483941449038287 = sqrt(2/pi)e^{-1/2} CORRECT; ratio to buggy = 2^{1/4} exactly.\n"
            "  GRreg 2^{9/4}=(2rt2)^{3/2} OK; GLm2reg 2^{19/4}=2^{9/4}·2^{5/2} OK; G0reg 2^{3/2} OK; TL0/mSreg/p2mag/d6 OK.\n"
            "  sympy: (t^2+rt2 t)^2+(rt2 t+1)^2=(t^2+rt2 t+1)^2 EXACT; H2 numerator 3x^2+6rt2 kx+4k^2; AM-GM (x-1)^2.\n")
print("logged")
