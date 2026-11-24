import pickle, sympy as sp
n = sp.Symbol('n')
Z = pickle.load(open('mu_zsum_K7.pkl','rb'))
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))
worst_deg = -99; all_poles = set(); budget_viol = []; sharp_viol = []; a1_viol = []
per_z = {}
for zkey, terms in sorted(Z.items()):
    md=-99; poles=set(); worstbudget=0
    for w, cs in terms:
        e = sp.together(sp.sympify(cs)); num, den = sp.fraction(e)
        dn = sp.Poly(sp.expand(num), n); dd = sp.Poly(sp.expand(den), n)
        dg = dn.degree()-dd.degree(); md=max(md,dg); worst_deg=max(worst_deg,dg)
        wl = sum(abs(x) for x in w)
        if wl > 7-wt(zkey): a1_viol.append((zkey,w,wl))
        if dd.degree()>0:
            for rt, mult in sp.roots(sp.expand(dd), n).items():
                rt = sp.nsimplify(rt); poles.add(rt); all_poles.add(rt)
                if 2*rt != int(2*rt): sharp_viol.append(('non-half-integer pole', zkey, w, rt))
                budget = wl + mult
                worstbudget = max(worstbudget, budget)
                if budget > 8 - wt(zkey): budget_viol.append((zkey, w, rt, mult))
                # sharp set {1/2, 0} u {-t/2: 1<=t<=6}
                ok_sharp = (rt == sp.Rational(1,2)) or (rt == 0) or (rt < 0 and -2*rt in range(1,7) and rt.denominator==2 or rt<0 and rt.is_integer)
                if not ok_sharp: sharp_viol.append((zkey, w, rt))
        else:
            worstbudget = max(worstbudget, wl)
            if wl > 8 - wt(zkey): budget_viol.append((zkey, w, 'constant-coeff'))
    per_z[zkey]=(md, sorted(poles,key=sp.N), worstbudget)
print('k=7 analysis per zeta-monomial (wt | maxdeg | max(|w|+ord) vs bound 8-wt | poles):')
for zkey,(md,poles,wb) in per_z.items():
    print(f'  z={zkey} wt={wt(zkey)}: deg={md} (proper<=-2: {md<=-2}), max(|w|+ord)={wb} (bound {8-wt(zkey)}), poles={poles}')
print()
print('ALL coefficients proper (deg<=-2):', worst_deg<=-2, '(worst', worst_deg,')')
print('all poles:', sorted(all_poles, key=sp.N))
print('(A1) |w| <= 7-wt(z) violations:', len(a1_viol))
print('budget |w|+ord <= 8-wt(z) violations (per word AND pole location):', len(budget_viol), budget_viol[:3])
print('sharp pole set {-t/2: t<=6} u {0,1/2} violations:', len(sharp_viol), sharp_viol[:3])
open('progress.log','a').write(f'step04c: k=7 Z-sum analysis: proper(deg<=-2)={worst_deg<=-2} worst={worst_deg}; poles={sorted(all_poles,key=sp.N)}; budget viol={len(budget_viol)}; sharp-set viol={len(sharp_viol)}; A1 viol={len(a1_viol)}\n')
