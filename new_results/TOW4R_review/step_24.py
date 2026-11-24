
import os, sys, pickle
import sympy as sp
work = "<project>/research/08_theory/TOW4R_review"
d, m = sp.symbols('d m'); i = sp.I
ef = pickle.load(open(work + "/tow4_copy/e_fits.pkl", 'rb'))
print("e_fits keys:", sorted(ef.keys()) if hasattr(ef,'keys') else type(ef))
def get(k):
    v = ef[k]
    return sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) if isinstance(v, str) else v
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out
if (7,3) in ef:
    v = get((7,3))
    print("e_fits[(7,3)] =", sp.expand(v))
    print("== -2i(d-5)/Pi3 ?", sp.expand(v - (-2*i*(d-5)/Pi(3))))==0 if True else None
    print("equal to stale fit:", sp.simplify(sp.expand(v + 2*i*(d-5)/Pi(3)))==0)
# now check displayed F(7,1),F(7,2),F(8,1..3) against pkl
raw = pickle.load(open(work + "/tow4_copy/forced_F78.pkl", 'rb'))
F = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) for k, v in raw.items()}
disp = {}
disp[(7,1)] = ( 16*i*d**7 - 416*i*d**6 + (80+4840*i)*d**5 - 160*d**4*m + (-1360-32520*i)*d**4 + 2720*d**3*m
   + (8980+135069*i)*d**3 + (-17960+240*i)*d**2*m + (-26960-342004*i)*d**2 - 240*i*d*m**2
   + (53920-1200*i)*d*m + (30300+479005*i)*d + 1200*i*m**2 - 60600*m - 278550*i ) / (480*Pi(3))
disp[(7,2)] = ( 4*i*d**4 - 68*i*d**3 + (6+455*i)*d**2 - 12*d*m + (-30-1375*i)*d + 60*m + 1500*i ) / (6*Pi(3))
disp[(8,3)] = ( 20*d**3 - 300*d**2 + (1549-10*i)*d + 20*i*m - 2592 ) / (4*Pi(4))
disp[(8,2)] = ( 16*d**6 - 432*d**5 + (5032-48*i)*d**4 + 96*i*d**3*m + (-32056+720*i)*d**3 - 1440*i*d**2*m
   + (116233-3684*i)*d**2 + (48+7368*i)*d*m + (-223115+6228*i)*d - 48*m**2 - 12456*i*m
   + 173022 ) / (48*Pi(4))
for k, v in disp.items():
    print(k, "displayed == pkl:", sp.expand(v - F[k]) == 0)
