import pickle, sympy as sp, mpmath as mp, sys
mp.mp.dps = 50
wd = '<project>/research/08_theory/K1f_hexagon_order5'
sys.path.insert(0, wd)
jtab = {k: sp.sympify(v) for k,v in pickle.load(open(f'{wd}/jtab5.pkl','rb')).items()}
x = sp.Symbol('x'); z2,z3,z5 = sp.symbols('z2 z3 z5')
exec(open(f'{wd}/step_31.py').read().split("if __name__")[0])  # atom_eval machinery
AT = {sp.Symbol(n): n for n in ATOMWORDS}
def evalJ(e, xx):
    zv = {z2: mp.zeta(2), z3: mp.zeta(3)}
    tot = mp.mpf(0)
    for term in sp.Add.make_args(e):
        co, fs = term.as_coeff_mul()
        val = mp.mpf(str(co))
        for f in fs:
            base, ex = (f.args if f.is_Pow else (f, 1)); ex = int(ex)
            if base == x: val *= xx**ex
            elif base in zv: val *= zv[base]**ex
            else: val *= atom_eval(ATOMWORDS and str(base), xx)**ex
        tot += val
    return tot
print("x->0 behaviour of non-local J entries (expect ~(-x)^c * const):")
for (n,a,b,c) in [(2,1,1,2),(0,0,0,5),(3,2,1,1),(4,1,0,3)]:
    e = jtab[(n,a,b,c)]
    v1 = evalJ(e, mp.mpf('0.001')); v2 = evalJ(e, mp.mpf('0.0001'))
    ratio = v1/v2
    print(f"  J({n};{a},{b},{c}): J(1e-3)={mp.nstr(v1,4)} J(1e-4)={mp.nstr(v2,4)} ratio={mp.nstr(ratio,4)} (expect ~10^{c})")
