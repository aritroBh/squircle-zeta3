import os, json, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
open('progress.log','w').write("DEC1: decagon P'''(1) closed form — progress log\n")

# ---- exact geometry ----
r5 = sp.sqrt(5)
phi = (1+r5)/2
s36 = sp.sqrt(10-2*r5)/4   # sin(pi/5)
R = 1/(2*phi**2*s36)       # circumradius
# vertices V0 (kink2, angle 162deg), V1 (kink3, angle -162deg): edge x = -h
h = 1/(2*phi)              # apothem = R cos(pi/10)?? check numeric below
halfL = R/(2*phi)          # half edge length = R sin(pi/10)
V0 = sp.Matrix([-h,  halfL])
V1 = sp.Matrix([-h, -halfL])
A0 = 1/(2*phi**2); A1 = 1/(2*phi**2)
B1 = 1/(2*phi**2); kap1 = phi-1
B4 = 1/(2*phi);    kap4 = 1/phi-1
B5 = 1/(2*phi)
print("h =", sp.simplify(h), " halfL =", sp.simplify(halfL))
L = 2*halfL
print("L =", sp.simplify(L), " =", sp.nsimplify(L))
print("L numeric:", mp.nstr(mp.mpf(sp.N(L,30)),25))
# exact check of term values along edge
t = sp.symbols('t')
y = V0 + t*(V1-V0)
ns = [sp.Matrix([sp.cos(sp.pi*i/5), sp.sin(sp.pi*i/5)]) for i in range(1,6)]
# trig exact: use cos/sin values
cos36 = phi/2; sin36 = s36; cos72 = (r5-1)/4; sin72 = sp.sqrt(10+2*r5)/4
nvals = [(cos36,sin36),(cos72,sin72),(-cos72,sin72),(-cos36,sin36),(-1,0)]
terms = []
for (nx,ny) in nvals:
    v = sp.simplify(nx*y[0]+ny*y[1])
    terms.append(v)
    print("term:", v)
# Expect: [B1(1+kap1 t), A0 t, A1(1-t), B4(1+kap4 t), B5] up to signs (abs)
print("check a1 = B1(1+kap1 t):", sp.simplify(terms[0] - B1*(1+kap1*t)))
print("check a2 = A0 t:", sp.simplify(terms[1] - A0*t))
print("check a3 = -A1(1-t)? :", sp.simplify(terms[2] + A1*(1-t)))
print("check a4 = -B4(1+kap4 t)?:", sp.simplify(terms[3] + B4*(1+kap4*t)))
print("check a5 = -B5?:", sp.simplify(terms[4] + B5))
with open('progress.log','a') as f:
    f.write("step: exact geometry. R=1/(2 phi^2 sin36), h=1/(2phi), L=1/(2 phi^3 sin36)=R/phi.\n"
            " terms: a1=B1(1+kap1 t), kap1=phi-1; a2=A0 t; a3=A1(1-t); a4=B4(1+kap4 t), kap4=1/phi-1; a5=B5 const.\n"
            " A0=A1=B1=1/(2phi^2), B4=B5=1/(2phi). ratios {phi,1/phi,1}; x1=-kap1=-1/phi, x4=-kap4=1/phi^2.\n")