
import os, sys, pickle
import sympy as sp
work = "<project>/research/08_theory/TOW4R_review"
sys.path.insert(0, work + "/tow4_copy")
d, m = sp.symbols('d m')
i = sp.I
raw = pickle.load(open(work + "/tow4_copy/forced_F78.pkl", 'rb'))
F = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) for k, v in raw.items()}
print("keys:", sorted(F.keys()))
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out
# zigzag numbers
A = {0:1,1:1,2:1,3:2,4:5,5:16,6:61,7:272}
def e_top(delta):
    if delta % 2 == 0:
        s = delta//2
        return (-1)**s * A[s+1] / Pi(s)
    s = (delta-1)//2
    return (-1)**(s+1) * i * A[s+1] * (d-s-2) / Pi(s)
for dl in [6,7,8]:
    print(f"F({dl},{dl//2}) =", sp.expand(F[(dl,dl//2)]), "  e_top =", sp.expand(e_top(dl)),
          "  equal:", sp.expand(F[(dl,dl//2)] - e_top(dl))==0)
# displayed claim: F(7,3) = -2i(d-5)/Pi(3)
disp73 = -2*i*(d-5)/Pi(3)
print("displayed F(7,3)==pkl?", sp.expand(F[(7,3)]-disp73)==0, "| e_top(7)=5i(d-5)/Pi3 == pkl?", sp.expand(F[(7,3)]-e_top(7))==0)
