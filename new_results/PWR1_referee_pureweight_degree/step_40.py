import pwr1_engine as E, sympy as sp, json, time
E.NMAX = 10  # extend
t0 = time.time()
tot2, pure2, B2 = E.total_log(2)
print("engine n=2 through tau^10 in", round(time.time()-t0,1), "s")
# re-verify vs D16 json through tau^10
js = json.load(open('../D16_pure_weight/tau_series_2d.json'))
def to_mono_dict(expr):
    out = {}
    for t in sp.expand(expr).as_ordered_terms() if expr != 0 else []:
        coeff, zpart = t.as_coeff_Mul()
        key = []
        for f in sp.Mul.make_args(zpart):
            if f == 1: continue
            b, e = f.as_base_exp()
            if b in E.Z.values():
                zi = [k for k,v in E.Z.items() if v==b][0]
                key.extend([zi]*int(e))
            else: key.append(str(f))
        k2 = ' '.join(map(str, sorted(key)))
        out[k2] = out.get(k2, 0) + coeff
    return {k: sp.nsimplify(v) for k, v in out.items()}
ok = True
for k in range(0, 11):
    mine = to_mono_dict(pure2.coeff(k))
    theirs = {kk: sp.Rational(vv) for kk, vv in js['logP'].get(str(k), {}).items()}
    if mine != theirs: ok=False; print("MISMATCH at", k, mine, theirs)
print("logP exact vs D16 through tau^10:", ok)
import pickle
pickle.dump({k: str(pure2.coeff(k)) for k in range(11)}, open('w_exact_n2_tau10.pkl','wb'))
print("w9 =", pure2.coeff(9))
print("w10 =", pure2.coeff(10))
