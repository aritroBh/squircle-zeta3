
import pickle, json, sympy as sp, mpmath as mp, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps=110
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
geo = json.load(open(f"{WORK}/geometry.json"))
eng = json.load(open(f"{WORK}/../K1d_engine_fixed/engine_arcs_fixed.json")) if False else json.load(open('research/08_theory/K1d_engine_fixed/engine_arcs_fixed.json'))
l2_ = mp.log(2); z2v = mp.pi**2/6; z3v = mp.zeta(3)
Li3half = 7*z3v/8 - z2v*l2_/2 + l2_**3/6
A1,B,kap,q0,V02,L2,lA0,lA1,lB = sp.symbols('A1 B kap q0 V02 L2 lA0 lA1 lB')
def coeff_num(name, e):
    num,den = clean[name]
    subs = {A1:sp.Rational(e['A1']), B:sp.Rational(e['B']), kap:sp.Rational(e['kap']), L2:sp.Integer(e['L2']),
            q0:sp.Integer(e['q0']), V02:sp.Integer(e['V02'])}
    return sp.simplify(num.subs(subs)/den.subs(subs))
print("=== Prop.7 reduction: cz3 from general formula vs K1d engine ===")
for bkey in ["0","1"]:
    for e in geo[bkey]:
        r=e['r']
        cL3 = coeff_num('Li3(x)',e); cL3L = coeff_num('Li3(x/(x-1))',e); cz = coeff_num('zeta3',e)
        if r=='2':   z3L3, z3L3L = sp.Rational(-3,4), sp.Rational(7,8)
        elif r=='1/2': z3L3, z3L3L = sp.Rational(7,8), sp.Rational(-3,4)
        else: z3L3 = z3L3L = sp.Integer(0)
        cz3_mine = sp.simplify(cz + cL3*z3L3 + cL3L*z3L3L)
        k1d = eng[f"b{bkey}_e{e['edge']}"]
        target = sp.sympify(k1d['cz3'])
        print(f"b={bkey} e{e['edge']} r={r}: my cz3 = {cz3_mine}  | K1d cz3 = {target} | diff = {sp.simplify(cz3_mine-target)}  (N_e={k1d['Ne']})")
print()
print("=== Landen identity check (110 dps) ===")
for xv in ['-2','-1','1/2','2/3','-1/3','3/4']:
    x=mp.mpf(xv)
    lhs = mp.polylog(3,x)+mp.polylog(3,1-x)+mp.polylog(3,1-1/x)
    rhs = z3v + z2v*mp.log(1-x) - mp.log(x)*mp.log(1-x)**2/2 + mp.log(1-x)**3/6
    print(f"x={xv:>4}: |lhs-rhs| = {mp.nstr(abs(lhs-rhs),3)}")
print("=== inversion identity Li3(1/y) for y=x/(x-1) ===")
for xv in ['-2','2/3']:
    x=mp.mpf(xv); y=x/(x-1)
    lhs = mp.polylog(3,1/y)+mp.polylog(3,y)
    rhs = -mp.log(-y)**3/6 - (mp.pi**2/6)*mp.log(-y)   # Li3(1/y)+Li3(y) = -ln^3(-y)/6 - zeta2 ln(-y)? check sign conventions
    print(f"x={xv:>4}: Li3(1/y)+Li3(y) = {mp.nstr(lhs,10)}, rhs try = {mp.nstr(rhs,10)}, |diff|={mp.nstr(abs(lhs-rhs),3)}")
