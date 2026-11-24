
import pickle, json, sympy as sp, mpmath as mp, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps = 120
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
geo = json.load(open(f"{WORK}/geometry_fresh.json"))
L2s = sp.Symbol('L2')
pars_syms = sp.symbols('A1 B kap q0 V02 lA0 lA1 lB')
atom_syms = sp.symbols('l3a l3l l2l1 l1_3 l2a l1_2 l1 z3 z2l1 z2 one')
basis_order = ['Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','ln^3(1-x)','Li2(x)','ln^2(1-x)','ln(1-x)','zeta3','zeta2*ln(1-x)','zeta2','1']
Eexpr = sum(num/L2s**sp.Rational(5,2)*atom for (num,den),atom in zip([clean[n] for n in basis_order], atom_syms))
f = sp.lambdify(pars_syms + (L2s,) + atom_syms, Eexpr, 'mpmath')
def logpq(rat):
    rat = sp.Rational(rat)
    primes = {}
    for p in [2,3,5,13]:
        e = sp.factorint(rat.p).get(p,0) - sp.factorint(rat.q).get(p,0)
        if e: primes[p]=e
    assert sp.prod(sp.Integer(p)**e for p,e in primes.items()) == rat
    return sum(e*mp.log(p) for p,e in primes.items())
z2v = mp.pi**2/6; z3v = mp.zeta(3)
results={}
for e in geo:
    kapv = sp.Rational(e['kap']); rv = e['r']
    x = 1-mp.mpf(rv)
    l1 = -mp.log(1-x); l2a = mp.polylog(2,x); l3a = mp.polylog(3,x); l3l = mp.polylog(3, x/(x-1))
    la0 = logpq(sp.Rational(e['A0'])/sp.Rational(e['w0']))
    la1 = logpq(sp.Rational(e['A1'])/sp.Rational(e['w1']))
    lb  = logpq(sp.Rational(e['B'])/sp.Rational(e['wo']))
    val = f(mp.mpf(sp.Rational(e['A1'])), mp.mpf(sp.Rational(e['B'])), mp.mpf(kapv),
            mp.mpf(sp.Integer(e['q0'])), mp.mpf(sp.Integer(e['V02'])), la0, la1, lb,
            mp.mpf(sp.Integer(e['L2'])), l3a, l3l, l2a*l1, l1**3, l2a, l1**2, l1, z3v, z2v*l1, z2v, mp.mpf(1))
    results[e['edge']] = val
    print(f"edge {e['edge']} (r={rv}, x={mp.nstr(x,4)}): E_arc = {mp.nstr(val,40)}")
tot = 2*sum(results.values())
print("P'''(1) total = 2*(sum 3 arcs) =", mp.nstr(tot,60))
json.dump({str(k): mp.nstr(v,110) for k,v in results.items()}, open(f"{WORK}/fresh_formula.json","w"))
mp.mp.dps=60
print("atom values used: Li2(-2)=",mp.nstr(mp.polylog(2,-2),12)," Li3(-2)=",mp.nstr(mp.polylog(3,-2),12)," Li3(2/3)=",mp.nstr(mp.polylog(3,sp.mpf(2)/3),12))
