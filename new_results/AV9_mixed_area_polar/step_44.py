
import subprocess, threading
def run_pslq8():
    r = subprocess.run(['<project>/.venv/bin/python','pslq2.py','8','300'],capture_output=True,text=True,timeout=1500)
    print("PSLQ c8@300:", r.stdout, r.stderr[-150:])
th = threading.Thread(target=run_pslq8, daemon=True)
th.start()

clean = r'''
import sympy as sp, pickle
b = {eval(k): sp.sympify(v) for k,v in pickle.load(open('scratchpad/circle_bivariate.pkl','rb')).items()}
# normalize log(2^k) -> k ln2 : sympy writes log(bignum); detect power-of-two
def norm(e):
    e = sp.expand(e)
    out = 0
    for term in sp.Add.make_args(e):
        term2 = term
        for f in sp.Mul.make_args(term):
            if f.func is sp.log and f.args[0].is_Integer:
                n = int(f.args[0])
                if n > 0 and (n & (n-1)) == 0:
                    k = n.bit_length()-1
                    term2 = term2.subs(f, k*sp.log(2))
        out += term2
    return sp.expand(out)
z3, z5 = sp.zeta(3), sp.zeta(5)
print("clean [zeta(3)] parts of b_jk:")
for key in sorted(b):
    e = norm(b[key])
    c3 = sp.expand(e).coeff(z3,1)
    if c3: print(f"  b_{key}: [zeta3] = {c3}   full = {e if key in [(0,3),(1,2),(2,2)] else '...'}")
print("\nfull clean key coefficients:")
for key in [(0,1),(1,1),(0,3),(1,2),(2,2),(0,4),(1,3),(0,5)]:
    print(f"b_{key} =", norm(b[key]))
'''
open('clean.py','w').write(clean)
r = subprocess.run(['<project>/.venv/bin/python','clean.py'],capture_output=True,text=True,timeout=600)
print(r.stdout); print(r.stderr[-300:])
th.join(timeout=1550)
