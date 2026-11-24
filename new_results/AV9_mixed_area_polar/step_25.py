
prt = r'''
import sympy as sp, pickle
b = {eval(k): sp.sympify(v) for k,v in pickle.load(open('scratchpad/circle_bivariate.pkl','rb')).items()}
z3, z5, L2, pi = sp.zeta(3), sp.zeta(5), sp.log(2), sp.pi
print("total-order | coefficient | zeta3-part | zeta5-part | has ln2?")
for n in range(1,6):
    for j in range(n+1):
        k = n-j
        e = sp.expand(b[(j,k)])
        c3 = sp.simplify(e.coeff(z3,1)) if e.has(z3) else 0
        c5 = sp.simplify(e.coeff(z5,1)) if e.has(z5) else 0
        print(f"b_({j},{k})  z3: {c3}   z5: {c5}   ln2: {e.has(L2)}")
# zeta3 polynomial at total order 3 along s2 = lam*s1
lam = sp.symbols('lam')
P3 = sum(b[(j,3-j)].coeff(z3,1)*lam**(3-j) for j in range(4))/pi
print("\n[zeta3]/pi along s2=lam*s at order 3:", sp.factor(P3))
print("roots:", sp.solve(P3, lam))
# order 4 zeta3 polynomial
P4 = sum(b[(j,4-j)].coeff(z3,1)*lam**(4-j) for j in range(5))/pi
print("[zeta3]/pi order 4:", sp.factor(P4))
print("self-dual lam=-1 order4 z3:", sp.simplify(P4.subs(lam,-1)))
# order 5 zeta3, zeta5
P5 = sum(b[(j,5-j)].coeff(z3,1)*lam**(5-j) for j in range(6))/pi
print("[zeta3]/pi order 5:", sp.factor(P5))
Q5 = sum(b[(j,5-j)].coeff(z5,1)*lam**(5-j) for j in range(6))/pi
print("[zeta5]/pi order 5:", sp.factor(Q5))
print("self-dual lam=-1: P5(-1)=", sp.simplify(P5.subs(lam,-1)), " Q5(-1)=", sp.simplify(Q5.subs(lam,-1)))
'''
open('prt.py','w').write(prt)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','prt.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-500:])
