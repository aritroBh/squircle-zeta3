import sympy as sp
p, y = sp.symbols('p y', positive=True)
g = -sp.log(1-y**p)/p
H = sp.exp(g)
code_lines = []
for k in [1,2,3]:
    Hk = sp.expand(sp.diff(H,p,k).subs(p,2))
    s = sp.pycode(Hk)
    code_lines.append(f"def _H{k}raw(y):\n    return {s}\n")
code = "\n".join(code_lines)
open('Hgen.py','w').write(code)
print(code[:1500])
