
import mpmath as mp, json, pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps = 110
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
geo = json.load(open(f"{WORK}/geometry_fresh.json"))
direct = json.load(open(f"{WORK}/fresh_direct.json"))
L2s = sp.Symbol('L2')
pars_syms = sp.symbols('A1 B kap q0 V02 lA0 lA1 lB')
atom_syms = sp.symbols('l3a l3l l2l1 l1_3 l2a l1_2 l1 z3 z2l1 z2 one')
basis_order = ['Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','ln^3(1-x)','Li2(x)','ln^2(1-x)','ln(1-x)','zeta3','zeta2*ln(1-x)','zeta2','1']
Eexpr = sum(num/L2s**sp.Rational(5,2)*atom for (num,den),atom in zip([clean[n] for n in basis_order], atom_syms))
f = sp.lambdify(pars_syms + (L2s,) + atom_syms, Eexpr, 'mpmath')
def logpq(rat):
    rat = sp.Rational(rat); out=mp.mpf(0)
    for p in [2,3]:
        e = sp.factorint(rat.p).get(p,0) - sp.factorint(rat.q).get(p,0)
        out += e*mp.log(p)
    assert sp.prod(sp.Integer(p)**(sp.factorint(rat.p).get(p,0) - sp.factorint(rat.q).get(p,0)) for p in [2,3])==rat
    return out
z2v=mp.pi**2/6; z3v=mp.zeta(3)
print(f"{'edge':>4} {'r':>4} {'residual |formula-direct|':>28}")
worst=mp.mpf(0); resids={}
for e in geo:
    rv=e['r']; x=1-mp.mpf(rv)
    l1=-mp.log(1-x); l2a=mp.polylog(2,x); l3a=mp.polylog(3,x); l3l=mp.polylog(3,x/(x-1))
    la0=logpq(sp.Rational(e['A0'])/sp.Rational(e['w0'])); la1=logpq(sp.Rational(e['A1'])/sp.Rational(e['w1'])); lb=logpq(sp.Rational(e['B'])/sp.Rational(e['wo']))
    val=f(mp.mpf(sp.Rational(e['A1'])),mp.mpf(sp.Rational(e['B'])),mp.mpf(sp.Rational(e['kap'])),
          mp.mpf(sp.Integer(e['q0'])),mp.mpf(sp.Integer(e['V02'])),la0,la1,lb,mp.mpf(sp.Integer(e['L2'])),
          l3a,l3l,l2a*l1,l1**3,l2a,l1**2,l1,z3v,z2v*l1,z2v,mp.mpf(1))
    d=abs(val-mp.mpf(direct[str(e['edge'])]))
    resids[e['edge']]=d; worst=max(worst,d)
    print(f"{e['edge']:>4} {rv:>4} {mp.nstr(d,3):>28}   formula={mp.nstr(val,25)}")
tf = 2*sum(mp.mpf(resids[k]) for k in resids)  # not meaningful; recompute totals
totf = 2*sum([f(mp.mpf(sp.Rational(e['A1'])),mp.mpf(sp.Rational(e['B'])),mp.mpf(sp.Rational(e['kap'])),mp.mpf(sp.Integer(e['q0'])),mp.mpf(sp.Integer(e['V02'])),logpq(sp.Rational(e['A0'])/sp.Rational(e['w0'])),logpq(sp.Rational(e['A1'])/sp.Rational(e['w1'])),logpq(sp.Rational(e['B'])/sp.Rational(e['wo'])),mp.mpf(sp.Integer(e['L2'])),mp.polylog(3,1-mp.mpf(e['r'])),mp.polylog(3,(1-mp.mpf(e['r']))/(-mp.mpf(e['r']))),mp.polylog(2,1-mp.mpf(e['r']))*(-mp.log(mp.mpf(e['r']))),(-mp.log(mp.mpf(e['r'])))**3,mp.polylog(2,1-mp.mpf(e['r'])),mp.log(mp.mpf(e['r']))**2,-mp.log(mp.mpf(e['r'])),z3v,-z2v*mp.log(mp.mpf(e['r'])),z2v,mp.mpf(1)) for e in geo])
print("worst per-arc residual:", mp.nstr(worst,3))
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_21-24 (Task 3): fresh hexagon V=(1,0),(3,3),(0,1)+negatives, w=(3,5,3), k=(1/18,1/18,1/18), ratios r={{3,1/3,1}} (x=-2, 2/3, 0). Formula vs INDEPENDENT direct d^3 jet-integration (own polar jet code + theta(t) change of variables, 100 dps): worst per-arc residual {mp.nstr(worst,3)}. Verifier vjet polar agrees too (exactly formula values).\n")
print("done")
