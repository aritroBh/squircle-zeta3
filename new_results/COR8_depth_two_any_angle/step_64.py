
import pickle, json, sys, re
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 120
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
terms = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/identity_terms.pkl','rb'))
def Lin(n, g):
    return mp.re(mp.polylog(n, g))
def RHSnDT(cc):
    tot = DT(cc)
    for kind, data, cf in terms:
        a = sum(mp.mpf(v.numerator)/v.denominator * cc**E for E, v in cf.items())
        if kind == 'Li4':
            tot += a * Lin(4, mp.mpf(sp.N(data.subs(c, cc), 40)))
        else:
            pv = mp.mpf(1)
            for k, arg in data:
                if k == 'log':
                    pv *= mp.log(mp.mpf(sp.N(arg.subs(c, cc), 40)))
                else:
                    pv *= Lin(int(k[2]), mp.mpf(sp.N(arg.subs(c, cc), 40)))
            tot += a * pv
    return tot
cc = mp.mpf('0.5')
Dval = Omega(cc) - RHSnDT(cc)
print("residual:", mp.nstr(Dval, 40))
# level-2 weight-4 constants basis
L2 = mp.log(2)
basis = [mp.zeta(4), mp.zeta(3)*L2, mp.zeta(2)*L2**2, L2**4, mp.polylog(4, mp.mpf('0.5')),
         mp.zeta(3), mp.zeta(2)*L2, L2**3, mp.polylog(3, mp.mpf('0.5')),
         mp.zeta(2), L2**2, mp.polylog(2, mp.mpf('0.5')), L2, mp.mpf(1)]
rel = mp.pslq([Dval] + basis, tol=mp.mpf('1e-110'), maxcoeff=10**8)
print("PSLQ:", rel)
