
# step_04.py -- COR6X referee: rational points 1/3, 2/5, 3/7, 1/2 (exact symbols, delta22)
#            + Task 3: Lie class of Omega(1/3) mod products
import sys; sys.path.insert(0, '.')
from cor6x_core import *
import sympy as sp, json, re

path = "<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"
data = json.load(open(path))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[1], sp.Rational(c)) for s, c in data if len(parse_term(s)[1]) == 4]
print("weight-4 terms:", len(TERMS))
assert len(TERMS) == 56

def symbol_Omega(c0):
    """exact symbol of Omega(c0) = sum coeff * H_w(c0), letters expanded in prime basis.
       (factors c0^e are weight-0 and drop out of the symbol)."""
    tot = {}
    for w, coeff in TERMS:
        tot = sym_add(tot, sym_scale(sym_H(w, c0), coeff))
    return tot

def primes_in(S):
    ps = set()
    for mono in S:
        for t in mono:
            ps |= set(expand_letter(t))
    return sorted(ps, key=str)

def delta22_primes(S):
    return delta22(sym_expanded(S))

for c0 in [sp.Rational(1,3), sp.Rational(2,5), sp.Rational(3,7), sp.Rational(1,2)]:
    S = symbol_Omega(c0)
    ps = primes_in(S)
    D = delta22_primes(S)
    print(f"\n=== c0 = {c0} ===")
    print("  symbol terms:", len(S), " prime letters:", ps, " letter-span dim:", len(ps))
    print("  delta22:", "ZERO" if not D else f"NONZERO, {len(D)} coordinates")
    if D:
        for k, v in sorted(D.items(), key=str):
            print("     ", k, v)

# ---- Task 3: Lie class of Omega(1/3) modulo products ----
print("\n=== Task 3: Omega(1/3) vs -(11/81)Li4(-8) - 2Li4(-3) + (20/3)Li4(2/3) mod products ===")
S_om = sym_expanded(symbol_Omega(sp.Rational(1,3)))
def sym_Li4(z):
    return sym_scale(sym_G((0,0,0,1), z), -1)
S_rhs = sym_expanded(sym_add(
    sym_scale(sym_Li4(-8), sp.Rational(-11,81)),
    sym_scale(sym_Li4(-3), sp.Rational(-2)),
    sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(20,3))))
diff = sym_sub(S_om, S_rhs)
print("  symbol difference has", len(diff), "terms")

# product subspace in weight 4 over the primes appearing
primes = sorted(set(t for mono in sym_add(S_om, S_rhs) for t in mono), key=str)
print("  primes:", primes)
import itertools
basis_monos = list(itertools.product(primes, repeat=4))
idx = {m: i for i, m in enumerate(basis_monos)}
rows = []
def add_prod_rows(u, v):
    for m1 in itertools.product(primes, repeat=u):
        for m2 in itertools.product(primes, repeat=v):
            r = [0]*len(basis_monos)
            for w in shuffle_words(m1, m2):
                r[idx[w]] += 1
            rows.append(r)
add_prod_rows(1,3); add_prod_rows(2,2); add_prod_rows(3,1)
P = sp.Matrix(rows)
print("  product subspace dim:", P.rank())
dv = [diff.get(m, 0) for m in basis_monos]
M2 = P.T
aug = M2.row_join(sp.Matrix(dv))
print("  difference in product subspace:", M2.rank() == aug.rank())

# ---- bonus: paper's claim for Omega(1/2): -21 Li4(-2) + 13/2 Li4(1/3) + 13/2 Li4(2/3) ----
print("\n=== bonus check: Omega(1/2) vs -21Li4(-2)+(13/2)Li4(1/3)+(13/2)Li4(2/3) mod products ===")
S_om2 = sym_expanded(symbol_Omega(sp.Rational(1,2)))
S_rhs2 = sym_expanded(sym_add(
    sym_scale(sym_Li4(-2), sp.Rational(-21)),
    sym_scale(sym_Li4(sp.Rational(1,3)), sp.Rational(13,2)),
    sym_scale(sym_Li4(sp.Rational(2,3)), sp.Rational(13,2))))
diff2 = sym_sub(S_om2, S_rhs2)
primes2 = sorted(set(t for mono in sym_add(S_om2, S_rhs2) for t in mono), key=str)
basis2 = list(itertools.product(primes2, repeat=4)); idx2 = {m:i for i,m in enumerate(basis2)}
rows2 = []
for m1 in itertools.product(primes2, repeat=1):
    for m2 in itertools.product(primes2, repeat=3):
        r = [0]*len(basis2)
        for w in shuffle_words(m1, m2): r[idx2[w]] += 1
        rows2.append(r)
for m1 in itertools.product(primes2, repeat=2):
    for m2 in itertools.product(primes2, repeat=2):
        r = [0]*len(basis2)
        for w in shuffle_words(m1, m2): r[idx2[w]] += 1
        rows2.append(r)
P2 = sp.Matrix(rows2)
dv2 = [diff2.get(m, 0) for m in basis2]
print("  primes:", primes2, " product subspace dim:", P2.rank(),
      " difference in product subspace:", P2.T.rank() == P2.T.row_join(sp.Matrix(dv2)).rank())
