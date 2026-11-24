import pwr1_engine as E, sympy as sp, json
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
        k2 = ' '.join(map(str, sorted(key, key=lambda x: int(x) if str(x).isdigit() else 999)))
        out[k2] = out.get(k2, 0) + coeff
    return {k: sp.nsimplify(v) for k, v in out.items()}

tot2, pure2, B2 = E.total_log(2)
P = E.ser_exp(pure2).scale(4)
ok = True
for name, ser in [('P/sqrt2', P), ('logP', pure2), ('B2', B2)]:
    jsdict = js[{'P/sqrt2':'P_over_sqrt2','logP':'logP','B2':'B2'}[name]]
    good = True
    for k in range(0, 9):
        mine = to_mono_dict(ser.coeff(k))
        theirs = {kk: sp.Rational(vv) for kk, vv in jsdict.get(str(k), {}).items()}
        if mine != theirs:
            good = ok = False
            print(f"  MISMATCH {name} tau^{k}: mine={mine}\n     theirs={theirs}")
    print(f"{name}: EXACT match vs D16 json through tau^8: {good}")
print("ALL:", ok)
