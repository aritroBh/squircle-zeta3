
import os, re, pickle
import sympy as sp
work = "<project>/research/08_theory/TOW4R_review"
tow4 = "<project>/research/08_theory/TOW4_tower_depth8"
d, m = sp.symbols('d m'); i = sp.I
N7 = (64*i*d**10 - 2048*i*d**9 + (672+29808*i)*d**8 - 1344*d**7*m + (-17472-259952*i)*d**7
 + 34944*d**6*m + (196560+1503964*i)*d**6 + (-393120+6720*i)*d**5*m + (-1241520-6034952*i)*d**5
 - 6720*i*d**4*m**2 + (2483040-114240*i)*d**4*m + (4758978+17057277*i)*d**4
 + 114240*i*d**3*m**2 + (-9514596+734160*i)*d**3*m + (-11077248-33707943*i)*d**3
 + (-10080-734160*i)*d**2*m**2 + (22137696-2113440*i)*d**2*m + (14484330+44731447*i)*d**2
 + 6720*d*m**3 + (50400+2113440*i)*d*m**2 + (-28968660+2293200*i)*d*m + (-8162700-35621565*i)*d
 - 33600*m**3 - 2293200*i*m**2 + 16325400*m + 12303900*i)
a7_disp = N7/(40320*(2*d-7)*(2*d-5)*(2*d-3))
a7pkl = sp.sympify(pickle.load(open(work+"/tow4_copy/a7_closed.pkl",'rb')), locals={'d':d,'m':m,'I':sp.I})
print("a7 displayed == pkl:", sp.expand(a7_disp - a7pkl) == 0)
# also check a7 == F(7,0)
raw = pickle.load(open(work + "/tow4_copy/forced_F78.pkl", 'rb'))
F = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) for k, v in raw.items()}
print("a7 displayed == F(7,0):", sp.expand(a7_disp - F[(7,0)]) == 0)
# a8 leading terms + denominator
a8pkl = sp.sympify(pickle.load(open(work+"/tow4_copy/a8_closed.pkl",'rb')), locals={'d':d,'m':m,'I':sp.I})
print("a8 displayed denom check: 645120*(2d-9)(2d-7)(2d-5)(2d-3)")
print("a8 pkl denom:", sp.denom(sp.together(a8pkl)).factor())
print("a8 == F(8,0):", sp.expand(a8pkl - F[(8,0)]) == 0)
num8 = sp.Poly(sp.numer(sp.together(a8pkl)), d, m)
print("a8 numer total degree in d (m=0):", sp.Poly(sp.numer(sp.together(a8pkl)).subs(m,0), d).degree(), "nterms:", len(sp.Poly(sp.numer(sp.together(a8pkl)), d, m).terms()))
lead = sp.Poly(sp.numer(sp.together(a8pkl)).subs(m, 0), d)
print("a8 leading terms:", [(t[0][0], t[1]) for t in lead.terms()[:3]])
