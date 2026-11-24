step_46 = r'''
import pickle, sympy as sp, json
d = pickle.load(open('Phi_syms.pkl','rb'))
PhiA, PhiL = d['PhiA'], d['PhiL']
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
x = sp.Symbol('x')
# evaluate at x = -kap (i.e. atom point 1-r) and factor
PhiA_k = sp.factor(PhiA.subs(x, -kap))
PhiL_k = sp.factor(PhiL.subs(x, -kap))
print("Phi_A(x=-kap) factored:")
print(PhiA_k)
print()
print("Phi_L(x=-kap) factored:")
print(PhiL_k)
print()
# collect by denominator structure; express as (B/kap^4)*poly?
for name, P in [('PhiA',PhiA_k),('PhiL',PhiL_k)]:
    num, den = sp.fraction(P)
    print(name, "denominator:", sp.factor(den))
    print(name, "numerator terms:", len(sp.expand(num).args) if sp.expand(num).is_Add else 1)
pickle.dump({'PhiA_k':PhiA_k,'PhiL_k':PhiL_k}, open('Phi_kappa.pkl','wb'))

# also with constraint A0 = A1 - B*kap substituted (single-ratio edge identity)
PhiA_c = sp.factor(PhiA_k.subs(A0, A1 - B*kap))
PhiL_c = sp.factor(PhiL_k.subs(A0, A1 - B*kap))
print("\nwith A0 = A1 - B*kap:")
print("Phi_A:", PhiA_c)
print()
print("Phi_L:", PhiL_c)
'''
open('step_46.py','w').write(step_46)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_46.py'],
                   capture_output=True, text=True, timeout=1800)
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
