
import pickle, json, sympy as sp
d = pickle.load(open('Phi_syms.pkl','rb'))
PhiA, PhiL = d['PhiA'], d['PhiL']   # NO factor 6 yet
A0,A1,B,kap,L2,q0,V02 = sp.symbols('A0 A1 B kap L2 q0 V02')
x = sp.Symbol('x')
geo = json.load(open('geometry.json'))

# candidate closed forms at x = -kap, using constraints A0=A1-B*kap, A1+B=1 (substitute A0 and B=1-A1)
G = (q0**2 - L2*V02)/L2**sp.Rational(5,2)
PA = 3*A0*kap*q0 - L2*(2*A0 + (kap+1)*A1)
QL = kap*(-A1**2 + A1*B*(2*kap-1) + 2*B**2*kap)*q0 + L2*(A0*B + (kap+1)*A1**2)
candA = B*G*PA
candL = 3*B*G*QL

# verify symbolically: substitute constraints into the DIFFERENCE
diffA = sp.expand(PhiA.subs(x, -kap) - candA)
diffL = sp.expand(PhiL.subs(x, -kap) - candL)
con = {A0: A1 - B*kap, B: 1 - A1}
print("PhiA(x=-kap) - candA, under constraints:", sp.expand(diffA.subs(con)))
print("PhiL(x=-kap) - candL, under constraints:", sp.expand(diffL.subs(con)))

# Now the FULL law with engine factor 6: N_e = PhA*[z3 of A(1-r)] + PhL*[z3 of Li3(1-r)]
z3cont = {'2':   {'A': sp.Rational(-5,8), 'L': sp.Rational(-3,4)},
          '1/2': {'A': sp.Rational(1),    'L': sp.Rational(7,8)}}
print()
print(f"{'edge':>12} {'PhA (6x)':>22} {'PhL (6x)':>22} {'N_e by law':>18} {'N_e step_42':>16}")
eng = json.load(open('engine_arcs_fixed.json'))
for bkey in ['0','1']:
    for e in geo[bkey]:
        if sp.Rational(e['kap']) == 0:
            print(f"b={bkey} e{e['edge']} r=1: N_e = 0 (kappa=0: W=log 1, atoms A(0)=Li3(0)=0)")
            continue
        p = {A0: sp.Rational(e['A0']), A1: sp.Rational(e['A1']), B: sp.Rational(e['B']),
             kap: sp.Rational(e['kap']), L2: sp.Integer(e['L2']), q0: sp.Integer(e['q0']),
             V02: sp.Integer(e['V02'])}
        PhA = sp.simplify(6*candA.subs(p)); PhL = sp.simplify(6*candL.subs(p))
        zc = z3cont[e['r']]
        Ne_law = sp.simplify(PhA*zc['A'] + PhL*zc['L'])
        Ne42 = sp.sympify(eng[f"b{bkey}_e{e['edge']}"]['Ne'])
        print(f"b={bkey} e{e['edge']} r={e['r']:>3}: {str(PhA):>22} {str(PhL):>22} {str(Ne_law):>18} {str(Ne42):>16}  match: {sp.simplify(Ne_law-Ne42)==0}")
