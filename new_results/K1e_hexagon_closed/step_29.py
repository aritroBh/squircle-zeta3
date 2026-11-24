
import pickle, json, sympy as sp, mpmath as mp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
mp.mp.dps = 130
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
geo = json.load(open(f"{WORK}/geometry.json"))
verif = {"0": json.load(open(f"{WORK}/t3_b0_150.json")), "1": json.load(open(f"{WORK}/t3_b1_110.json"))}
L2s = sp.Symbol('L2')
basis_order = ['Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','ln^3(1-x)','Li2(x)','ln^2(1-x)','ln(1-x)','zeta3','zeta2*ln(1-x)','zeta2','1']
# atom values at x=-1 (r=2) and x=1/2 (r=1/2), 130 digits
l2_ = mp.log(2); z2v = mp.pi**2/6; z3v = mp.zeta(3)
Li3half = 7*z3v/8 - z2v*l2_/2 + l2_**3/6
ATOM = {
 '2':   {'l1':-l2_, 'l2a':-z2v/2, 'l3a':-3*z3v/4, 'l3l':Li3half},
 '1/2': {'l1': l2_, 'l2a': z2v/2 - l2_**2/2, 'l3a':Li3half, 'l3l':-3*z3v/4},
}
def logpq(rat):
    rat = sp.Rational(rat)
    a2 = sp.factorint(rat.p).get(2,0) - sp.factorint(rat.q).get(2,0)
    a3 = sp.factorint(rat.p).get(3,0) - sp.factorint(rat.q).get(3,0)
    assert sp.Rational(2)**a2 * sp.Rational(3)**a3 == rat
    return a2*mp.log(2) + a3*mp.log(3)
def eval_arc(e):
    A1v = sp.Rational(e['A1']); Bv = sp.Rational(e['B']); kapv = sp.Rational(e['kap'])
    L2v = sp.Integer(e['L2']); q0v = sp.Integer(e['q0']); V02v = sp.Integer(e['V02'])
    A0v = sp.Rational(e['A0'])
    w0,w1,wo = sp.Rational(e['w0']), sp.Rational(e['w1']), sp.Rational(e['wo'])
    subs = {sp.Symbol('A1'):A1v, sp.Symbol('B'):Bv, sp.Symbol('kap'):kapv, L2s:L2v,
            sp.Symbol('q0'):q0v, sp.Symbol('V02'):V02v,
            sp.Symbol('lA0'):0, sp.Symbol('lA1'):0, sp.Symbol('lB'):0}
    # numeric coefficients: do lA's numerically since they are logs
    la0, la1, lb = logpq(A0v/w0), logpq(A1v/w1), logpq(Bv/wo)
    tot = mp.mpf(0)
    for name in basis_order:
        num, den = clean[name]
        cn = sp.N(num.subs(subs), 50)
        cd = sp.N(den.subs(subs), 50)
        # substitute logs numerically: recompute with lA numeric
        tot0 = mp.mpf(str(cn))/mp.mpf(str(cd))
        tot += tot0  # placeholder; we redo below properly
    return None
# better: build one big lambda
pars_syms = sp.symbols('A1 B kap q0 V02 lA0 lA1 lB')
Eexpr = sum(num/L2s**sp.Rational(5,2)*atom for (num,den),atom in
            zip([clean[n] for n in basis_order],
                [sp.Symbol(a) for a in ['l3a','l3l','l2l1','l1_3','l2a','l1_2','l1','z3','z2l1','z2','one']]))
atom_syms = sp.symbols('l3a l3l l2l1 l1_3 l2a l1_2 l1 z3 z2l1 z2 one')
f = sp.lambdify(pars_syms + (L2s,) + atom_syms, Eexpr, 'mpmath')
rows=[]
for bkey in ["0","1"]:
    arcs = verif[bkey]['arcs']
    for e in geo[bkey]:
        A1v=float(sp.Rational(e['A1'])); Bv=float(sp.Rational(e['B'])); kapv=float(sp.Rational(e['kap']))
        L2v=float(sp.Integer(e['L2'])); q0v=float(sp.Integer(e['q0'])); V02v=float(sp.Integer(e['V02']))
        A0v=sp.Rational(e['A0']); w0,w1,wo=sp.Rational(e['w0']),sp.Rational(e['w1']),sp.Rational(e['wo'])
        la0,la1,lb = logpq(A0v/w0), logpq(A1v/w1), logpq(Bv/wo)
        r = e['r']
        if kapv == 0:
            at = dict(l1=mp.mpf(0),l2a=mp.mpf(0),l3a=mp.mpf(0),l3l=mp.mpf(0))
        else:
            at = ATOM[r]
        val = f(mp.mpf(A1v),mp.mpf(Bv),mp.mpf(kapv),mp.mpf(q0v),mp.mpf(V02v),la0,la1,lb,mp.mpf(L2v),
                at['l3a'],at['l3l'],at['l2a']*at['l1'],at['l1']**3,at['l2a'],at['l1']**2,at['l1'],
                z3v,z2v*at['l1'],z2v,mp.mpf(1))
        P3 = mp.mpf(arcs[e['edge']-1]['P3_arc'])
        diff = abs(val-P3)
        rows.append((bkey,e['edge'],r,val,P3,diff))
        print(f"b={bkey} edge {e['edge']} (r={r}): E_arc={mp.nstr(val,30)}")
        print(f"   P3_arc={mp.nstr(P3,30)}  |diff|={mp.nstr(diff,3)}")
worst = max(rw[5] for rw in rows)
print("worst |diff| over 6 arcs:", mp.nstr(worst,3))
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_11: numerical check of general 11-basis closed form vs verifier P3_arc on ALL 6 arcs (b=0 150d, b=1 110d) at 130 dps: worst |diff| = {mp.nstr(worst,3)}.\n")
