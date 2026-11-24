import pickle, sys, time
sys.path.insert(0, ".")
from qs3_symbol import q3, ZERO, ONE, add, scale, shuffle, const_pow
import dod5d_grewrite as GR
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
GT = pickle.load(open("gsym_tab.pkl","rb"))
CA_VEC   = (0, -2, 0, 0); CC0_VEC  = (-1, 0, 0, 0); CC1_VEC  = (0, -2, 1, 0)
from collections import defaultdict
t0 = time.time()
byconst = defaultdict(lambda: {})
blocks = defaultdict(lambda: {})
def block_of(le):
    nz = tuple(i for i in range(4) if le[2+i] > 0)
    return 'local' if not nz else (f'single-{nz[0]}' if len(nz)==1 else (f'double-{nz[0]}{nz[1]}' if len(nz)==2 else (f'triple-{nz[0]}{nz[1]}{nz[2]}' if len(nz)==3 else 'quad')))
for (n, le), coef in E5X.items():
    if sum(le) != 5: continue
    vwt = le[0]+le[1]+le[2]+le[3]+le[4]+le[5]
    rw = rewrite[(n, le[0], le[1], le[2], le[3], le[4], le[5])]
    top = {}
    for w, cw in rw.items():
        if len(w) != vwt: continue
        top = add(top, scale(GT[w], cw))
    top = scale(top, coef)
    cpat = (le[6], le[7], le[8])
    byconst[cpat] = add(byconst[cpat], top)
    blocks[(cpat, block_of(le))] = add(blocks[(cpat, block_of(le))], top)
t1 = time.time()
def const_sym(e6, e7, e8):
    out = {(): ONE}
    if e6: out = shuffle(out, const_pow(CA_VEC, e6))
    if e7: out = shuffle(out, const_pow(CC0_VEC, e7))
    if e8: out = shuffle(out, const_pow(CC1_VEC, e8))
    return out
total = {}; block_syms = defaultdict(lambda: {})
for cpat, top in byconst.items():
    s = shuffle(top, const_sym(*cpat)) if any(cpat) else top
    total = add(total, s)
for (cpat, blk), top in blocks.items():
    s = shuffle(top, const_sym(*cpat)) if any(cpat) else top
    block_syms[blk] = add(block_syms[blk], s)
total = scale(total, q3(1440))
block_syms = {k: scale(v, q3(1440)) for k, v in block_syms.items()}
t2 = time.time()
from qs3_symbol import MONE
D = pickle.load(open("dod5d_symbols.pkl","rb"))
ok_total = add(total, scale(D['total'], MONE)) == {}
ok_blocks = all(add(block_syms[k], scale(D['blocks'][k], MONE)) == {} for k in D['blocks']) and set(block_syms) == set(D['blocks'])
print("grouped assembly: %.1fs + %.1fs; total==stored: %s; blocks==stored: %s" % (t1-t0, t2-t1, ok_total, ok_blocks))
