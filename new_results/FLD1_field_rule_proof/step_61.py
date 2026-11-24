WD='<project>/research/08_theory/FLD1_field_rule_proof'
import pickle, sympy as sp
n=sp.Symbol('n')
Z = pickle.load(open(WD+'/mu_zsum_K5.pkl','rb'))
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))
def poleorder(cs):
    e = sp.together(sp.sympify(cs)); num,den = sp.fraction(e)
    den = sp.expand(den)
    if sp.Poly(den,n).degree()==0: return 0
    # max multiplicity of any root
    rts = sp.roots(den, n)
    return max(rts.values()) if rts else 0
print("Full weight budget: max over terms of (|w| + maxpoleorder(rho_w)) <= k+1-wt(z)?")
allok=True
for k in range(1,6):
    for zkey, terms in sorted(Z[k].items()):
        worst = 0
        for w,c in terms:
            po = poleorder(c)
            val = sum(abs(x) for x in w) + po
            worst = max(worst, val)
        bound = k+1-wt(zkey)
        ok = worst <= bound
        allok = allok and ok
        print(f'  k={k} z={zkey}: max(|w|+poleorder) = {worst}, bound k+1-wt(z) = {bound}, ok={ok}')
print('ALL weight budgets ok:', allok)