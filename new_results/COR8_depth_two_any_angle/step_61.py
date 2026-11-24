
import pickle, sys
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 50
from fractions import Fraction as Fr
import sympy as sp
c = sp.Symbol('c')
terms = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/identity_terms.pkl','rb'))
def Lin(n, g):
    return mp.re(mp.polylog(n, g))
cc = mp.mpf('0.5')
rows = []
for kind, data, cf in terms:
    a = sum(mp.mpf(v.numerator)/v.denominator * cc**E for E, v in cf.items())
    if kind == 'Li4':
        g = data.subs(c, cc)
        v = Lin(4, mp.mpf(sp.N(g, 40)))
        contrib = a*v
        rows.append((abs(contrib), f"a={mp.nstr(a,8)} Li4({sp.sstr(data)}) = Li4({mp.nstr(mp.mpf(sp.N(g,20)),10)}) = {mp.nstr(v,15)} -> {mp.nstr(contrib,10)}"))
    else:
        pv = mp.mpf(1); desc = []
        for k, arg in data:
            av = arg.subs(c, cc)
            if k == 'log':
                x = mp.log(mp.mpf(sp.N(av, 40)))
            else:
                x = Lin(int(k[2]), mp.mpf(sp.N(av, 40)))
            pv *= x
            desc.append(f"{k}({mp.nstr(mp.mpf(sp.N(av,20)),8)})={mp.nstr(x,10)}")
        contrib = a*pv
        rows.append((abs(contrib), f"a={mp.nstr(a,8)} {' '.join(desc)} -> {mp.nstr(contrib,10)}"))
for mag, s in sorted(rows, reverse=True)[:14]:
    print(s)
