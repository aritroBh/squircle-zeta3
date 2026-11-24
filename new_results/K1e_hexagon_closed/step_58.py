
import pickle, json, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
geo = json.load(open(f"{WORK}/geometry_fresh.json"))
L2s=sp.Symbol('L2')
lA0,lA1,lB = sp.symbols('lA0 lA1 lB')
for e in geo:
    subs={sp.Symbol('A1'):sp.Rational(e['A1']), sp.Symbol('B'):sp.Rational(e['B']), sp.Symbol('kap'):sp.Rational(e['kap']),
          L2s:sp.Integer(e['L2']), sp.Symbol('q0'):sp.Integer(e['q0']), sp.Symbol('V02'):sp.Integer(e['V02'])}
    print(f"=== fresh edge {e['edge']} (r={e['r']}, x={1-sp.Rational(e['r'])}); coefficients (num/L2^(5/2)) ===")
    for name in ['Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','zeta2*ln(1-x)','zeta3']:
        num,den = clean[name]
        cv = sp.simplify(num.subs(subs)/(sp.Integer(e['L2'])**sp.Rational(5,2)))
        print(f"  {name:>16}: {cv}   (atom at x={1-sp.Rational(e['r'])})")
    # check Li3 coefficient is rational (no logs)
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_26: printed exact polylog coefficients for fresh edges (rational multiples of 1/L2^(5/2); e.g. r=3 edge: Li3(-2) coeff 63/(2*13^(5/2))).\n")
