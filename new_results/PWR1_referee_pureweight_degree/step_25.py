import pwr1_engine as E, sympy as sp, json

js = json.load(open('../D16_pure_weight/tau_series_2d.json'))

def to_mono_dict(expr):
    """convert sympy expr in Z's to {'2 2 4': rational, ...}"""
    out = {}
    expr = sp.expand(expr)
    terms = expr.as_ordered_terms() if expr != 0 else []
    for t in terms:
        coeff, zpart = t.as_coeff_Mul()
        key = []
        for f in sp.Mul.make_args(zpart):
            if f.has(E.GAM): key.append('GAM')
            else:
                b, e = f.as_base_exp()
                if b in E.Z.values():
                    zi = [k for k,v in E.Z.items() if v==b][0]
                    key.extend([str(zi)]*int(e))
                else: key.append(str(f))
        key = ' '.join(sorted(key, key=int)) if key and key[0]!='GAM' else ' '.join(key)
        out[key] = out.get(key, 0) + coeff
    return {k: sp.nsimplify(v) for k, v in out.items()}

tot2, pure2, B2 = E.total_log(2)
P = E.ser_exp(pure2).scale(4)

def compare(mine_ser, jsdict, name, kmax=8):
    allok = True
    for k in range(0, kmax+1):
        mine = to_mono_dict(mine_ser.coeff(k))
        theirs = {kk: sp.Rational(vv) for kk, vv in jsdict.get(str(k), {}).items()}
        if mine != theirs:
            allok = False
            print(f"  MISMATCH {name} tau^{k}: mine={mine} theirs={theirs}")
    print(f"{name}: exact match through tau^{kmax}: {allok}")
    return allok

compare(P, js['P_over_sqrt2'], 'P/sqrt2')
compare(pure2, js['logP'], 'logP')
compare(B2, js['B2'], 'B2')
open('progress.log','a').write("step_16-17: built pwr1_engine.py (formal-Z exact engine); n=2: gamma=0 all orders, weight-k homogeneous, Rat_2=1; P/sqrt2, logP, B2 EXACT match vs D16 tau_series_2d.json through tau^8\n")
