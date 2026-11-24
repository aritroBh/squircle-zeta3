
import sys, itertools, pickle
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
Z4 = (Fr(0),)*4
def v4(j, val=1):
    v = list(Z4); v[j] = Fr(val); return tuple(v)
L2, cv, mv, pv = v4(0), v4(1), v4(2), v4(3)
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/symbols_2kept.pkl','rb'))
R, polefree = D['R'], D['polefree']
# ansatz
li4 = {s: li_sym(4, rv, omv) for s, rv, omv in polefree}
logs = {'c': cv, '1-c': mv, '1+c': pv}
prods = {}
log2terms = {}
for s, rv, omv in polefree:
    S3 = li_sym(3, rv, omv); S2 = li_sym(2, rv, omv)
    log2terms[f"l2*Li3({s})"] = sshuffle(log_sym(L2), S3)
    log2terms[f"l2^2*Li2({s})"] = sshuffle(sshuffle(log_sym(L2), log_sym(L2)), S2)
    for hn, hv in logs.items():
        prods[f"Li3({s})*log({hn})"] = sshuffle(S3, log_sym(hv))
        log2terms[f"l2*Li2({s})*log({hn})"] = sshuffle(log_sym(L2), sshuffle(S2, log_sym(hv)))
        for h2n, h2v in logs.items():
            prods[f"Li2({s})*log({hn})*log({h2n})"] = sshuffle(S2, sshuffle(log_sym(hv), log_sym(h2v)))
    prods[f"Li2({s})**2"] = sshuffle(S2, S2)
for i, (s1, rv1, omv1) in enumerate(polefree):
    for s2, rv2, omv2 in polefree[i+1:]:
        prods[f"Li2({s1})*Li2({s2})"] = sshuffle(li_sym(2, rv1, omv1), li_sym(2, rv2, omv2))
for combo in itertools.combinations_with_replacement(list(logs), 4):
    S = log_sym(logs[combo[0]])
    for h in combo[1:]:
        S = sshuffle(S, log_sym(logs[h]))
    prods["*".join(f"log({h})" for h in combo)] = S
for combo in itertools.combinations_with_replacement(list(logs), 3):
    S = log_sym(L2)
    for h in combo:
        S = sshuffle(S, log_sym(logs[h]))
    log2terms["l2*" + "*".join(f"log({h})" for h in combo)] = S
for combo in itertools.combinations_with_replacement(list(logs), 2):
    S = sshuffle(log_sym(L2), log_sym(L2))
    for h in combo:
        S = sshuffle(S, log_sym(logs[h]))
    log2terms["l2^2*" + "*".join(f"log({h})" for h in combo)] = S
for hn, hv in logs.items():
    S = sshuffle(sshuffle(log_sym(L2), log_sym(L2)), log_sym(L2))
    log2terms[f"l2^3*log({hn})"] = sshuffle(S, log_sym(hv))
print("li4:", len(li4), " prods:", len(prods), " log2-terms:", len(log2terms))
# to 256-dim vectors
bwords = list(itertools.product(range(4), repeat=4))
bidx = {w: i for i, w in enumerate(bwords)}
def expand_letter(l):
    return [(j, l[j]) for j in range(4) if l[j]]
def to256(S, e=0):
    v = [Fr(0)]*256
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
def dropz(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
li4v = {k: to256(dropz(v)) for k, v in li4.items()}
prodv = {k: to256(dropz(v)) for k, v in prods.items()}
log2v = {k: to256(dropz(v)) for k, v in log2terms.items()}
Rv = {E: to256(R, E) for E in (-1, 0, 1)}
print("vectors built")
pickle.dump({'li4v': li4v, 'prodv': prodv, 'log2v': log2v, 'Rv': Rv}, open('solve256.pkl','wb'))
import numpy as np
p = (1<<31)-1
def rank_modp(vecs):
    basis = []
    for v in vecs:
        v = np.array([x % p for x in v], dtype=np.int64)
        for b, pw in basis:
            if v[pw]: v = (v - int(v[pw])*b) % p
        nz = np.nonzero(v)[0]
        if len(nz):
            pw = int(nz[0]); basis.append(((v*pow(int(v[pw]), p-2, p)) % p, pw))
    return len(basis)
allv = list(li4v.values()) + list(prodv.values()) + list(log2v.values())
ra = rank_modp(allv)
print("rank all:", ra, " R_E in span:", all(rank_modp(allv + [Rv[E]]) == ra for E in (-1,0,1)))
print("rank products:", rank_modp(list(prodv.values())), " products+log2:", rank_modp(list(prodv.values())+list(log2v.values())), " Li4 alone:", rank_modp(list(li4v.values())))
