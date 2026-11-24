import os, json, sympy as sp, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
r5 = sp.sqrt(5); phi = (1+r5)/2
t = sp.symbols('t')
s36sq = (5-r5)/8            # sin^2(pi/5)
s36 = sp.sqrt(s36sq)        # sin(pi/5)  (radical, only enters overall L)
R = 1/(2*phi**2*s36)
h = 1/(2*phi)
halfL2 = 1/(2*phi**6*(5-r5))   # halfL^2 in Q(sqrt5)
Lexpr = 2*sp.sqrt(halfL2)
print("L^2 =", sp.simplify(4*halfL2), "; L =", sp.simplify(Lexpr))
mp.mp.dps=40
print("L num:", mp.nstr(mp.mpf(sp.N(Lexpr,40)),30))
# term values (closed forms, verified numerically at 50 dps earlier)
A0=A1=B1=1/(2*phi**2); B4=B5=1/(2*phi); kap1=phi-1; kap4=1/phi-1
a = [A0*t, A1*(1-t), B1*(1+kap1*t), B4*(1+kap4*t), B5*(1+0*t)]
print("sum a =", sp.simplify(sum(a)))
# C_5 from E1 edge law: Z_e = 6 a b (h^2/L)(3(a u + b v)-2), a=b=1/(2phi^2), u=v=1/2
al = 1/(2*phi**2); be = al; u=sp.Rational(1,2); v=u
Ze = 6*al*be*(h**2/Lexpr)*(3*(al*u+be*v)-2)
C5 = sp.simplify(10*Ze)
print("Z_e =", sp.simplify(Ze))
print("C_5 =", C5, "=", sp.nsimplify(C5))
C5n = mp.mpf(sp.N(C5,50))
print("C_5 numeric:", mp.nstr(C5n,40))
z3c = json.load(open('research/08_theory/D15b_polygon_law/z3coeff_m2_8.json'))
print("D15b z3coeff keys:", list(z3c.keys()))
print("D15b C_5 (120d):", z3c['5'] if '5' in z3c else z3c)
# full 200-digit P'''(1)
v=json.load(open('research/08_theory/D15b_polygon_law/vals_m2_8_200.json'))
print("P3 200d:", v['5'][3])
print("P 200d:", v['5'][0])
diff = C5n - mp.mpf(z3c['5'][0] if isinstance(z3c['5'],list) else z3c['5'])
print("|C_5 - D15b| =", mp.nstr(abs(diff),6))
with open('progress.log','a') as f:
    f.write(f"step: Q(sqrt5) geometry set. L = {sp.simplify(Lexpr)}; overall radical only in L.\n"
            f" C_5 = 10 Z_e = {C5} ; |C_5 - D15b 120d value| = {mp.nstr(abs(diff),4)}\n")
pickle.dump({'A0':A0,'A1':A1,'B1':B1,'B4':B4,'B5':B5,'kap1':kap1,'kap4':kap4,'L':Lexpr,'h':h,'halfL2':halfL2,'R':R,'phi':phi,'C5':C5,'Ze':Ze}, open('geometry.pkl','wb'))
print("saved geometry.pkl")