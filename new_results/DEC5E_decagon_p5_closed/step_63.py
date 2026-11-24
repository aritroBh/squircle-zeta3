import sys, pickle, math, itertools
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib as D
from fractions import Fraction as F
QS5 = D.QS5
bs = pickle.load(open('../DEC5D_decagon_order5_depth/dec_block_symbols.pkl','rb'))
sub = bs['sub_sym']
# exact log part per block (assuming symmetric symbols)
import mpmath as mp
LOGPART = {}
for blk in ['local','single-1','single-4','mixed']:
    s = sub[blk]
    sym_ok = True
    for w in s:
        for p in set(itertools.permutations(w)):
            if p not in s or s[p] != s[w]: sym_ok = False; break
        if not sym_ok: break
    part = {}
    for a in range(6):
        w = tuple([0]*a + [2]*(5-a))
        c = s.get(w, QS5(0)) / QS5(math.factorial(a)*math.factorial(5-a))
        if c != QS5(0): part[a] = c
    LOGPART[blk] = part
    print(blk, "symmetric:", sym_ok, " log part:", {a: str(c) for a, c in part.items()})
pickle.dump({k: {a: (str(v.p), str(v.q)) for a, v in p.items()} for k, p in LOGPART.items()}, open('block_logparts.pkl','wb'))
# sanity: sum of block log parts = total log part
tot = {}
for blk, p in LOGPART.items():
    for a, c in p.items(): tot[a] = tot.get(a, QS5(0)) + c
print("sum of blocks:", {a: str(c) for a, c in tot.items()})
