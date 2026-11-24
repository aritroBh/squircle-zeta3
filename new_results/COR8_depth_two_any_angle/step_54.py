
import pickle, json, sys, re
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 80
from hpl_series import H_series, Li22
from fractions import Fraction as Fr
import sympy as sp
c = sp.Symbol('c')
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle/'
R3 = json.load(open(G1+'hpl_fit_k3_(0,0,0).json'))
om_terms = []
for name, coef in R3:
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    if len(w) == 4:
        om_terms.append((e, w, Fr(coef)))
def Omega(cc):
    return sum(mp.mpf(q.numerator)/q.denominator * cc**e * H_series(w, cc) for e, w, q in om_terms)
def DT(cc):
    return ((1-cc)/cc)*Li22(mp.mpf('0.5'), 1-cc) - ((1+cc)/cc)*Li22(mp.mpf('0.5'), 1+cc)
def Lin(n, g):
    return mp.re(mp.polylog(n, g))
terms = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/identity_terms.pkl','rb'))
def RHS(cc):
    tot = DT(cc)
    for kind, data, cf in terms:
        a = sum(mp.mpf(v.numerator)/v.denominator * cc**E for E, v in cf.items())
        if kind == 'Li4':
            tot += a * Lin(4, mp.mpf(sp.N(data.subs(c, cc), 30)))
        else:
            pv = mp.mpf(1)
            for k, arg in data:
                if k == 'log':
                    av = {'c': cc, '1-c': 1-cc, '1+c': 1+cc}[str(arg)]
                    pv *= mp.log(av)
                else:
                    pv *= Lin(int(k[2]), mp.mpf(sp.N(arg.subs(c, cc), 30)))
            tot += a * pv
    return tot
for cc in [mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('0.7'), mp.mpf('0.13')]:
    d = Omega(cc) - RHS(cc)
    print(f"c={cc}: |Omega-RHS| = {mp.nstr(abs(d),3)}   Omega={mp.nstr(Omega(cc),20)}")
