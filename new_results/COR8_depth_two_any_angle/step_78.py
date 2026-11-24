
import pickle, json, sys, re, time
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 55
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
# precompute HPL coefficient arrays once (float)
from hpl_series import hpl_coeffs
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
LI4EX = {nm: sp.sympify(nm) for (k_, nm) in alln if k_ == 'Li4'}
PRODEX = {nm: [(k, (sp.sympify(a) if k!='log' else a)) for k,a in v] for nm,v in PRODS.items()}
def Lin(n, g): return mp.re(mp.polylog(n, g))
used = [(j, kind, nm, {E: sols[E][j] for E in (-1,0,1) if sols[E][j]}) for j,(kind,nm) in enumerate(alln) if any(sols[E][j] for E in (-1,0,1))]
def resid(cc):
    tot = Omega(cc) - DT(cc)
    for j, kind, nm, cf in used:
        a = sum(mp.mpf(v.numerator)/v.denominator * cc**E for E, v in cf.items())
        if kind == 'Li4':
            tot -= a * Lin(4, mp.mpf(sp.N(LI4EX[nm].subs(c, cc), 40)))
        else:
            pv = mp.mpf(1)
            for k, arg in PRODEX[nm]:
                if k == 'log':
                    av = {'c': cc, '1-c': 1-cc, '1+c': 1+cc}[arg]
                    pv *= mp.log(av)
                else:
                    pv *= Lin(int(k[2]), mp.mpf(sp.N(arg.subs(c, cc), 40)))
            tot -= a * pv
    return tot
# invisible basis (rebuild same as before)
L2 = mp.log(2)
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
sa = json.load(open(C4+"sunit_args.json"))
args = [s for _,_,s in sa if s not in ('1','-1','2','1/2')]
argex = {s: sp.sympify(s) for s in args}
def lg(cc, i): return mp.log([cc, 1-cc, 1+cc][i])
def Li(n, ex, cc): return mp.re(mp.polylog(n, mp.mpf(sp.N(ex.subs(c,cc),40))))
basis = []
basis.append(('z4', lambda cc: mp.zeta(4)))
basis.append(('z3*l2', lambda cc: mp.zeta(3)*L2))
basis.append(('z2*l2^2', lambda cc: mp.zeta(2)*L2**2))
basis.append(('l2^4', lambda cc: L2**4))
for g, i in [('c',0),('1-c',1),('1+c',2)]:
    basis.append((f'z3*log({g})', lambda cc, i=i: mp.zeta(3)*lg(cc,i)))
    basis.append((f'z2*l2*log({g})', lambda cc, i=i: mp.zeta(2)*L2*lg(cc,i)))
    basis.append((f'l2^3*log({g})', lambda cc, i=i: L2**3*lg(cc,i)))
import itertools
for s in args:
    ex = argex[s]
    basis.append((f'z2*Li2({s})', lambda cc, ex=ex: mp.zeta(2)*Li(2, ex, cc)))
    basis.append((f'l2^2*Li2({s})', lambda cc, ex=ex: L2**2*Li(2, ex, cc)))
    basis.append((f'l2*Li3({s})', lambda cc, ex=ex: L2*Li(3, ex, cc)))
    for h, j in [('c',0),('1-c',1),('1+c',2)]:
        basis.append((f'l2*Li2({s})*log({h})', lambda cc, ex=ex, j=j: L2*Li(2, ex, cc)*lg(cc,j)))
for (h1,i1),(h2,i2) in itertools.combinations_with_replacement([('c',0),('1-c',1),('1+c',2)],2):
    basis.append((f'z2*log({h1})*log({h2})', lambda cc, i1=i1, i2=i2: mp.zeta(2)*lg(cc,i1)*lg(cc,i2)))
    basis.append((f'l2^2*log({h1})*log({h2})', lambda cc, i1=i1, i2=i2: L2**2*lg(cc,i1)*lg(cc,i2)))
for trip in itertools.combinations_with_replacement([('c',0),('1-c',1),('1+c',2)],3):
    idx = [t[1] for t in trip]
    nm = '*'.join(f'log({t[0]})' for t in trip)
    basis.append((f'l2*{nm}', lambda cc, idx=idx: L2*mp.prod([lg(cc,i) for i in idx])))
print(len(basis), "basis")
t0 = time.time()
pts = [mp.mpf('0.04') + i*mp.mpf('0.0032') for i in range(253)]
# evaluate residual and basis at points
import numpy as np
rv = []
for p_ in pts:
    rv.append(resid(p_))
print("residual evals done", time.time()-t0)
B = []
for nm, f in basis:
    B.append([f(p_) for p_ in pts])
print("basis evals done", time.time()-t0)
pickle.dump({'rv': [str(x) for x in rv], 'B': [[str(x) for x in col] for col in B],
             'pts': [str(p_) for p_ in pts]}, open('fit_data.pkl','wb'))
