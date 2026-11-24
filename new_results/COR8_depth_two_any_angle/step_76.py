
import pickle, json, sys, re, time
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 70
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
S = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/full_sols_v2.pkl','rb'))
sols, alln = S['sols'], S['alln']
def parse_prod(nm):
    parts = []
    i = 0
    while i < len(nm):
        m2 = re.match(r'(Li[23]|log)\(', nm[i:])
        kind = m2.group(1)
        j = i + len(kind) + 1
        depth = 1; k = j
        while depth:
            if nm[k] == '(': depth += 1
            elif nm[k] == ')': depth -= 1
            k += 1
        arg = nm[j:k-1]
        parts.append((kind, arg))
        if nm[k:k+3] == '**2':
            parts.append((kind, arg)); k += 3
        if k < len(nm):
            assert nm[k] == '*'
            k += 1
        i = k
    return parts
PRODS = {nm: parse_prod(nm) for (k_, nm) in alln if k_ == 'prod'}
def Lin(n, g):
    return mp.re(mp.polylog(n, g))
def term_val(kind, nm, cc):
    if kind == 'Li4':
        return Lin(4, mp.mpf(sp.N(sp.sympify(nm).subs(c, cc), 40)))
    pv = mp.mpf(1)
    for k, arg in PRODS[nm]:
        if k == 'log':
            av = {'c': cc, '1-c': 1-cc, '1+c': 1+cc}[arg]
            pv *= mp.log(av)
        else:
            pv *= Lin(int(k[2]), mp.mpf(sp.N(sp.sympify(arg).subs(c, cc), 40)))
    return pv
def resid(cc):
    tot = Omega(cc) - DT(cc)
    for j, (kind, nm) in enumerate(alln):
        cf = {E: sols[E][j] for E in (-1,0,1) if sols[E][j]}
        if not cf: continue
        a = sum(mp.mpf(v.numerator)/v.denominator * cc**E for E, v in cf.items())
        tot -= a * term_val(kind, nm, cc)
    return tot
for cs in ['0.5','0.3']:
    print(f"residual({cs}) =", mp.nstr(resid(mp.mpf(cs)), 30))
pickle.dump('ok', open('/tmp/r.pkl','wb'))
