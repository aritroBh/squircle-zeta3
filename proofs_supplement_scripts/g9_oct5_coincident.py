# Check of the split-entry rule for coincident letters (g9_oct5_symbol.py, function symI): the symbol of G(a,a;1) must be
# the symbol of G(a;1)^2/2, and G(a)G(b)=G(a,b)+G(b,a), for the letters a,b in {0,p,q} (1 is excluded: G(1;1) diverges).
import os, sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g9_oct5_symbol as S
def same(d1, d2):
    keys = set(d1) | set(d2)
    return all((d1.get(k, S.Q2(0)) - d2.get(k, S.Q2(0))).iszero() for k in keys)
ok = True
for a in 'pq':
    sa = S.symG_words((a,)); saa = S.symG_words((a, a))
    prod = S.shuffle_sym(sa, sa)
    half = {k: v * S.Q2(2) for k, v in saa.items()}
    r = same(prod, half); ok &= r
    print("G(%s,%s;1): Sym(G(%s)^2) = 2 Sym G(%s,%s):" % (a, a, a, a, a), r)
for a, b in [('p', 'q'), ('q', 'p'), ('0', 'p'), ('p', '0'), ('0', 'q')]:
    sa, sb = S.symG_words((a,)), S.symG_words((b,))
    lhs = S.shuffle_sym(sa, sb); rhs = dict(S.symG_words((a, b)))
    for k, v in S.symG_words((b, a)).items():
        S.sadd(rhs, k, v)
    r = same(lhs, rhs); ok &= r
    print("G(%s)G(%s) = G(%s,%s)+G(%s,%s):" % (a, b, a, b, b, a), r)
assert ok
print('pass')
