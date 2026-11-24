# A5: transcription check: the closed forms DISPLAYED in proof_depth8.md  vs  my independently derived a_7,a_8,F(7,s),F(8,s)
import sys, re, pickle
sys.dont_write_bytecode=True
W='<project>/research/08_theory/TOW4_tower_depth8'
sys.path.insert(0, W)
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application
import tow2lib as T2
d,m = T2.d, T2.m; I = sp.I
txt = open(W+'/proof_depth8.md').read()
def P(s):
    s = s.replace('^','**').replace('\n',' ')
    s = re.sub(r'(?<![A-Za-z_])i(?![A-Za-z_])', 'I', s)
    s = s.replace('Pi_3','PI3').replace('Pi_4','PI4')
    return parse_expr(s, transformations=standard_transformations+(implicit_multiplication_application,), local_dict={'d':d,'m':m,'I':I,'PI3':T2.Pi(3),'PI4':T2.Pi(4)})
mine = {k: sp.sympify(v, locals={'d':d,'m':m,'I':I}) for k,v in pickle.load(open('a_mine.pkl','rb')).items()}
# N_7
i0 = txt.index('N_7 ='); i1 = txt.index('a_8 = N_8')
blk = txt[i0+len('N_7 ='):i1]
blk = blk.replace('/ (40320','')   # safety
N7 = P(blk.split('a_8')[0].strip())
a7disp = N7/(40320*T2.Pi(3))
print('displayed a_7 == mine a_7 :', sp.cancel(sp.together(a7disp - mine[7])) == 0)
# N_8 leading terms
den8 = 645120*T2.Pi(4)
N8 = sp.cancel(sp.together(mine[8]*den8)); N8p = sp.Poly(sp.expand(N8), d, m)
lead = sp.Poly(sp.expand(N8), d).all_coeffs()[:3]
print('N_8 mine leading coeffs in d:', [sp.simplify(c) for c in lead])
print('   displayed: 256 d^12 - 10752 d^11 + (208384 - 3584 i) d^10  ;  match:', [sp.simplify(c) for c in lead] == [256, -10752, 208384-3584*I],
      '  terms:', len(N8p.terms()), '(claimed 57)', ' N_7 terms:', len(sp.Poly(sp.expand(sp.cancel(mine[7]*40320*T2.Pi(3))), d, m).terms()), '(claimed 41)')
# pickles a7_closed/a8_closed (strings)
for k,fn in ((7,'a7_closed.pkl'),(8,'a8_closed.pkl')):
    raw = pickle.load(open(W+'/'+fn,'rb'))
    ex = sp.sympify(raw, locals={'d':d,'m':m,'I':I})
    print(fn, 'equals mine:', sp.cancel(sp.together(ex - mine[k])) == 0)
# F displays
try:
    F = {k: sp.sympify(v, locals={'d':d,'m':m,'I':I}) for k,v in pickle.load(open('F_mine.pkl','rb')).items()}
except Exception as e:
    print('F_mine.pkl not ready', e); sys.exit()
def grab(label, nxt_label_regex, pi, den):
    j = txt.index(label)
    rest = txt[j+len(label):]
    # numerator is between the first '(' and ') /'
    a = rest.index('(')
    b = rest.index(') /')
    num = P(rest[a:b+1])
    return num/(den*pi)
disp = {}
disp[(7,1)] = grab('F(7,1) =', None, T2.Pi(3), 480)
disp[(7,2)] = grab('F(7,2) =', None, T2.Pi(3), 6)
disp[(8,1)] = grab('F(8,1) =', None, T2.Pi(4), 5760)
disp[(8,2)] = grab('F(8,2) =', None, T2.Pi(4), 48)
disp[(8,3)] = grab('F(8,3) =', None, T2.Pi(4), 4)
disp[(7,3)] = -2*I*(d-5)/T2.Pi(3)
disp[(8,4)] = 16/T2.Pi(4)
for k,v in sorted(disp.items()):
    print(f'displayed F{k} == mine:', sp.cancel(sp.together(v - F[k])) == 0)
