
import pickle, json, sys
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 80
from hpl_series import H_series, Li22
from fractions import Fraction as Fr
import sympy as sp

D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/full_sols_raw.pkl','rb'))
sols, alln = D['sols'], D['alln']
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle/'
R3 = json.load(open(G1+'hpl_fit_k3_(0,0,0).json'))
c = sp.Symbol('c')

# Omega: weight-4 HPL terms of R3
import re
om_terms = []
for name, coef in R3:
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    if len(w) == 4:
        om_terms.append((e, w, Fr(coef)))
print("Omega terms:", len(om_terms))

def Omega(cc):
    tot = mp.mpf(0)
    for e, w, q in om_terms:
        tot += mp.mpf(q.numerator)/q.denominator * cc**e * H_series(w, cc)
    return tot

def DT(cc):
    return ((1-cc)/cc)*Li22(mp.mpf('0.5'), 1-cc) - ((1+cc)/cc)*Li22(mp.mpf('0.5'), 1+cc)

def Lin(n, g):
    v = mp.polylog(n, g)
    return mp.re(v) if abs(mp.im(v)) > 0 else mp.re(v)

def prodval(nm, cc):
    gv = lambda s: sp.sympify(s).subs(c, cc)
    if '**2' in nm:
        base = nm[4:-3]  # strip Li2( and )**2
        return Lin(2, mp.mpf(sp.N(gv(base), 30)))**2
    parts = nm.split('*')
    tot = mp.mpf(1)
    for p in parts:
        if p.startswith('Li3('):
            tot *= Lin(3, mp.mpf(sp.N(gv(p[4:-1]), 30)))
        elif p.startswith('Li2('):
            tot *= Lin(2, mp.mpf(sp.N(gv(p[4:-1]), 30)))
        elif p.startswith('log('):
            arg = {'c': cc, '1-c': 1-cc, '1+c': 1+cc}[p[4:-1]]
            tot *= mp.log(arg)
        else:
            raise ValueError(p)
    return tot

def RHS(cc):
    tot = DT(cc)
    for j, (kind, nm) in enumerate(alln):
        cf = {E: sols[E][j] for E in (-1,0,1) if sols[E][j]}
        if not cf: continue
        a = sum(mp.mpf(v.numerator)/v.denominator * cc**E for E, v in cf.items())
        if kind == 'Li4':
            tot += a * Lin(4, mp.mpf(sp.N(sp.sympify(nm).subs(c, cc), 30)))
        else:
            tot += a * prodval(nm, cc)
    return tot

for cc in [mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('0.7')]:
    d = Omega(cc) - RHS(cc)
    print(f"c={cc}: Omega={mp.nstr(Omega(cc),25)}  RHS={mp.nstr(RHS(cc),25)}  |diff|={mp.nstr(abs(d),3)}")
