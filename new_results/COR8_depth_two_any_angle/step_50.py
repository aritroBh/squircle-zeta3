
import sys, json, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
import itertools
Z3 = (Fr(0),)*3
def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)
zero = v3(0,0)
def drop_zero(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
bwords = list(itertools.product(range(3), repeat=4))
bidx = {w:i for i,w in enumerate(bwords)}
def expand_letter(l):
    return [(j, l[j]) for j in range(3) if l[j]]
def to81(S, e=0):
    v = [Fr(0)]*81
    for w, cf in S.items():
        coeff = cf.get(e)
        if coeff is None: continue
        terms = [((), Fr(1))]
        for l in w:
            opts = expand_letter(l)
            terms = [(p+(j,), c0*m) for p, c0 in terms for j, m in opts]
        for pp, c0 in terms:
            v[bidx[pp]] += coeff*c0
    return v
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
li4v, prodv, Rv_full = D['li4v'], D['prodv'], D['Rv']  # Rv was Omega - DT(a)
# need S(Omega) at 81-level: Rv + DT(a)
S22m = drop_zero(li22_sym(zero, mv, zero, cv, pv))
S22p = drop_zero(li22_sym(zero, pv, zero, cv, mv))
S31m = drop_zero(li31_sym(zero, mv, zero, cv, pv))
S31p = drop_zero(li31_sym(zero, pv, zero, cv, mv))
def vadd(u, v, s=1): return [x + s*y for x, y in zip(u, v)]
def vscale(u, q): return [q*x for x in u]
# DT symbols at 81 level, per exponent
def dt_vec(Sm, Sp, cm, cp):
    out = {}
    for E in (-1,0,1):
        v = vscale(to81(Sm), cm[E])
        v = vadd(v, vscale(to81(Sp), cp[E]))
        out[E] = v
    return out
variants = {
 'a_Li22pair': dt_vec(S22m, S22p, {-1:Fr(1),0:Fr(-1),1:Fr(0)}, {-1:Fr(-1),0:Fr(-1),1:Fr(0)}),
 'b_Li31pair': dt_vec(S31m, S31p, {-1:Fr(-1),0:Fr(1),1:Fr(0)}, {-1:Fr(1),0:Fr(1),1:Fr(0)}),
 'c_mixed31m22p': dt_vec(S31m, S22p, {-1:Fr(-1),0:Fr(1),1:Fr(0)}, {-1:Fr(-1),0:Fr(-1),1:Fr(0)}),
}
# Omega at 81 level = Rv + DT(a)
OM = {E: vadd(Rv_full[E], variants['a_Li22pair'][E]) for E in (-1,0,1)}
li4n = list(li4v.keys()); prodn = list(prodv.keys())
allv = [li4v[k] for k in li4n] + [prodv[k] for k in prodn]
def solve(Rres):
    M = [[allv[j][i] for j in range(len(allv))] + [Rres[-1][i], Rres[0][i], Rres[1][i]] for i in range(81)]
    basis = []
    for r in M:
        for b, pw in basis:
            if r[pw]:
                f = r[pw]; r = [x - f*y for x, y in zip(r, b)]
        nz = next((j for j in range(len(allv)) if r[j]), None)
        if nz is None: continue
        f = r[nz]; basis.append(([x/f for x in r], nz))
    sols = {}
    for k, E in enumerate((-1,0,1)):
        rhs = len(allv) + k
        sol = [Fr(0)]*len(allv)
        for b, pw in reversed(basis):
            sol[pw] = b[rhs] - sum(b[j]*sol[j] for j in range(pw+1, len(allv)) if sol[j])
        sols[E] = sol
    return sols
for name, dtv in variants.items():
    Rres = {E: vadd(OM[E], dtv[E], -1) for E in (-1,0,1)}
    sols = solve(Rres)
    used = set()
    for E in (-1,0,1):
        for j in range(len(allv)):
            if sols[E][j]: used.add(j)
    nli4 = sum(1 for j in used if j < len(li4n))
    print(f"{name}: total terms {len(used)} ({nli4} Li4, {len(used)-nli4} products)")
pickle.dump(OM, open('omega81.pkl','wb'))
