src = '''"""dod5d_symbol_assembly.py -- assemble the full weight-5 symbol of P^(5)(1)/u for the 12-gon.
For each E5 monomial ((n, le), coef=(p,q)) of total weight 5:
   top G-words (length = variable log weight) -> g_symbol(0, w, 1)
   constant logs CA (slot6), CC0 (slot7), CC1 (slot8) shuffled in via their V-logvecs
Everything scaled by 1440 (= 12 edges * 5!). Block = pattern of nonzero W_i slots.
Output: dod5d_symbols.pkl with per-block and total weight-5 symbols (slot words over V={0,1,2[,3]}).
"""
import pickle, sys, time
sys.path.insert(0, ".")
from qs3_symbol import (q3, Q3V, ZERO, ONE, add, scale, shuffle, g_symbol, const_pow)
import dod5d_grewrite as GR

E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))

CA_VEC   = (0, -2, 0, 0)     # log alpha = log((2-r3)/2)
CC0_VEC  = (-1, 0, 0, 0)     # log c0_0 = log(2-r3)
CC1_VEC  = (0, -2, 1, 0)     # log c0_1 = log((2r3-3)/2)

def block_of(le):
    nz = tuple(i for i in range(4) if le[2+i] > 0)
    if len(nz) == 0: return 'local'
    if len(nz) == 1: return f'single-{nz[0]}'
    if len(nz) == 2: return f'double-{nz[0]}{nz[1]}'
    if len(nz) == 3: return f'triple-{nz[0]}{nz[1]}{nz[2]}'
    return 'quad'

SLOT = {GR.L0: 0, GR.L1: 1, GR.P0: 2, GR.P1: 3, GR.P2: 4, GR.P3: 5}

t0 = time.time()
# group monomials by const-log pattern to share the shuffle
from collections import defaultdict
byconst = defaultdict(lambda: {})   # (e6,e7,e8) -> accumulated top symbol (unscaled by monomial coef yet)
blocks = defaultdict(lambda: {})
gsym_cache = {}
def GSYM(w):
    if w not in gsym_cache:
        gsym_cache[w] = g_symbol(ZERO, w, ONE)
    return gsym_cache[w]

nmon = 0
for (n, le), coef in E5X.items():
    wtot = sum(le)
    if wtot != 5: continue
    nmon += 1
    vwt = le[0]+le[1]+le[2]+le[3]+le[4]+le[5]
    rw = rewrite[(n, le[0], le[1], le[2], le[3], le[4], le[5])]
    top = {}
    for w, cw in rw.items():
        if len(w) != vwt: continue
        top = add(top, scale(GSYM(w), cw))
    top = scale(top, coef)
    cpat = (le[6], le[7], le[8])
    byconst[cpat] = add(byconst[cpat], top)
    blocks[(cpat, block_of(le))] = add(blocks[(cpat, block_of(le))], top)
print("top symbols assembled for %d monomials, %d const patterns, %d g_symbols, %.1fs" %
      (nmon, len(byconst), len(gsym_cache), time.time()-t0))

def const_sym(e6, e7, e8):
    out = {(): ONE}
    if e6: out = shuffle(out, const_pow(CA_VEC, e6))
    if e7: out = shuffle(out, const_pow(CC0_VEC, e7))
    if e8: out = shuffle(out, const_pow(CC1_VEC, e8))
    return out

total = {}
block_syms = defaultdict(lambda: {})
for cpat, top in byconst.items():
    s = shuffle(top, const_sym(*cpat)) if any(cpat) else top
    total = add(total, s)
for (cpat, blk), top in blocks.items():
    s = shuffle(top, const_sym(*cpat)) if any(cpat) else top
    block_syms[blk] = add(block_syms[blk], s)
# scale by 1440
total = scale(total, q3(1440))
block_syms = {k: scale(v, q3(1440)) for k, v in block_syms.items()}
print("total symbol entries:", len(total), " blocks:", {k: len(v) for k, v in sorted(block_syms.items())})
pickle.dump({'total': total, 'blocks': dict(block_syms)}, open("dod5d_symbols.pkl","wb"))
print("saved dod5d_symbols.pkl  (%.1fs)" % (time.time()-t0))
'''
with open("dod5d_symbol_assembly.py","w") as f: f.write(src)
import subprocess
r = subprocess.run(["<project>/.venv/bin/python", "-u", "dod5d_symbol_assembly.py"],
                   capture_output=True, text=True, timeout=1700)
print(r.stdout)
print(r.stderr[-1500:])
