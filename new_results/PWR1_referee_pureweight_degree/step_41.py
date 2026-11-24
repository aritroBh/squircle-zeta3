src = open('pwr1_engine.py').read()
src = src.replace('NMAX = 8\n', 'NMAX = 10\n')
open('pwr1_engine.py','w').write(src)
import importlib, pwr1_engine as E
importlib.reload(E)
import sympy as sp, json, time
t0=time.time()
tot2, pure2, B2 = E.total_log(2)
print("engine n=2 through tau^10 in", round(time.time()-t0,1),"s")
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
ok = all(to_mono_dict(pure2.coeff(k)) == {kk: sp.Rational(vv) for kk, vv in js['logP'].get(str(k), {}).items()} for k in range(11))
print("logP exact vs D16 through tau^10:", ok)
w_ok = all(E.weight_check(pure2.coeff(k), k) for k in range(11))
gam_ok = all(sp.expand(tot2.coeff(k)).coeff(E.GAM,1)==0 for k in range(11))
print("weight-homogeneous all k<=10:", w_ok, "; gamma-free:", gam_ok)
json.dump({str(k): str(pure2.coeff(k)) for k in range(11)}, open('w_exact_n2_tau10.json','w'))
print("w9 =", pure2.coeff(9))
