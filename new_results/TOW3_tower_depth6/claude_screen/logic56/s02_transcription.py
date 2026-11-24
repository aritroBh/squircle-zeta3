# s02: (a) tau-table vs definition; (b) printed closed forms in proof_depth6.md vs tow3lib; (c) pickles vs lib
import sys, re, time, pickle
sys.dont_write_bytecode = True
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
import tow2lib as T2, tow3lib as T3
d, m = T3.d, T3.m
I = sp.I
t0=time.time()
# (a) tau constants from the DL1/DL2 Taylor series definition (independent of tow2lib._tau_tab)
v = sp.symbols('v')
a = 1+sp.sqrt(2); b = sp.sqrt(2)-1
def taus(DL, nmax):
    out=[]
    w = sp.symbols('w')
    f = DL.subs(v, I+w)
    ser = sp.series(f, w, 0, nmax+3).removeO()
    for n in range(nmax+1):
        out.append(sp.nsimplify(sp.simplify(sp.expand(ser.coeff(w, n+1)))))
    return out
DL1 = 1/(v-b)+1/(v+a)-1/v
DL2 = 1/(v-a)+1/(v+b)-1/v
t1 = taus(DL1, 7); t2 = taus(DL2, 7)
ok=True
for n in range(8):
    e1 = sp.simplify(sp.radsimp(t1[n]) - T2.tau1f(n)); e2 = sp.simplify(sp.radsimp(t2[n]) - T2.tau2f(n))
    ok &= (e1==0 and e2==0)
    # conj relation conj(tau1_n) = (-1)^n tau2_n
    ok &= sp.simplify(sp.conjugate(T2.tau1f(n)) - (-1)**n*T2.tau2f(n))==0
print('(a) tau table n<=7 equals Taylor coefficients of DL1/DL2 at v=i (own series), and conj(tau1_n)=(-1)^n tau2_n:', ok, round(time.time()-t0,1),'s', flush=True)

# (b) printed forms in proof_depth6.md
md = open('<project>/research/08_theory/TOW3_tower_depth6/proof_depth6.md').read()
tr = standard_transformations + (implicit_multiplication_application, convert_xor)
loc = {'d': d, 'm': m, 'I': I}
def P(s):
    s = re.sub(r'\s+', ' ', s)
    s = re.sub(r'\bi\b', 'I', s)
    return parse_expr(s, local_dict=loc, transformations=tr)
def grab(start, end):
    i0 = md.index(start); i1 = md.index(end, i0)
    return md[i0+len(start):i1]
Pi = T3.Pi
a5_txt = grab('a_5 = (', ')\n        / (480')
a6_txt = grab('a_6 = (', ')\n        / (5760')
a5 = P(a5_txt)/(480*(2*d-5)*(2*d-3))
a6 = P(a6_txt)/(5760*(2*d-7)*(2*d-5)*(2*d-3))
print('(b1) printed a_5 == tow3lib.arho(5):', sp.simplify(a5 - T3.arho(5))==0, flush=True)
print('(b2) printed a_6 == tow3lib.arho(6):', sp.simplify(a6 - T3.arho(6))==0, round(time.time()-t0,1), flush=True)
pk = {k: sp.sympify(v, locals={'d':d,'m':m,'I':I}) for k,v in pickle.load(open('arho56_simplified.pkl','rb')).items()}
print('(b3) arho56_simplified.pkl == lib arho:', sp.simplify(pk[5]-T3.arho(5))==0, sp.simplify(pk[6]-T3.arho(6))==0, flush=True)
e51 = P(grab('e(5,1) = (', ')/(12 Pi_2)'))/(12*Pi(2))
e61 = P(grab('e(6,1) = (', ')/(96 Pi_3)'))/(96*Pi(3))
e62 = P(grab('e(6,2) = (', ')/(2 Pi_3)'))/(2*Pi(3))
e41 = P(grab('e(4,1) = (', ')/Pi_2'))/Pi(2)
print('(b4) printed e(4,1),e(5,1),e(6,1),e(6,2) == lib:', [sp.simplify(x-y)==0 for x,y in [(e41,T3.e41),(e51,T3.e51),(e61,T3.e61),(e62,T3.e62)]], flush=True)
print('    printed e(5,2) = -2i(d-4)/Pi_2, e(6,3) = -5/Pi_3, e(4,2)=2/Pi_2  vs lib e():',
      sp.simplify(T3.e(5,2) - (-2*I*(d-4)/Pi(2)))==0, sp.simplify(T3.e(6,3) - (-5/Pi(3)))==0, sp.simplify(T3.e(4,2)-2/Pi(2))==0)
# (c) pickles
ef = pickle.load(open('e_fits.pkl','rb'))
def fit_expr(rec):
    # monomials (i,j) -> d^i m^j with (re,im) pair coefficients, over Pi_b
    num = 0
    for (ex, ey), (re_, im_) in zip(rec['mons'], rec['sol']):
        num += (sp.Rational(re_.numerator, re_.denominator) + I*sp.Rational(im_.numerator, im_.denominator))*d**ex*m**ey
    return num/Pi(rec['b'])
for key, lib in [((4,1),T3.e41),((5,1),T3.e51),((6,2),T3.e62)]:
    print('(c) e_fits.pkl', key, '== lib:', sp.simplify(fit_expr(ef[key])-lib)==0, flush=True)
print('    e61_closed.pkl str:', pickle.load(open('e61_closed.pkl','rb'))[:200])
print('done', round(time.time()-t0,1))
