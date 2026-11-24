import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, add, scale, shuffle, gsym, const_pow, obstruction
from fractions import Fraction as F

SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(SRC+"parsed_mons.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))

t0 = time.time()
# symbol of each needed integral's TOP weight part: sum_{w, |w|=a+b+c+d} gamma_w S(G(w;1))
symcache = {}
def Gsym(w):
    if w not in symcache: symcache[w] = gsym(w)
    return symcache[w]

def top_symbol(a,b,c,d,n):
    k = a+b+c+d
    out = {}
    for w, gam in rewrite[(n,a,b,c,d)].items():
        if len(w) != k: continue
        out = add(out, scale(Gsym(w), gam))
    return out

def qs2(p,q): return Q(F(p.numerator,p.denominator), F(q.numerator,q.denominator))

block_sym = {'local': {}, 'single': {}, 'mixed': {}}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    for n,(p,q) in cu.items():
        coef = qs2(p,q)
        sI = top_symbol(a,b,c,d,n)
        # shuffle with (log2)^e (log(1+sqrt2))^f symbols
        s = sI
        if e: s = shuffle(s, const_pow(0, e))
        if f: s = shuffle(s, const_pow(1, f))
        block_sym[cls] = add(block_sym[cls], scale(s, coef))
block_sym = {k: scale(v, Q(960)) for k,v in block_sym.items()}
print("block symbols computed in %.1fs" % (time.time()-t0))
for k,v in block_sym.items():
    print(k, "symbol words:", len(v))
# save
def ser(s): return {''.join(map(str,w)): (str(c.p), str(c.q)) for w,c in s.items()}
pickle.dump({k: ser(v) for k,v in block_sym.items()}, open("block_symbols.pkl","wb"))

# local block symbol should be exactly 0 (letters {0,1}, differences ±1)
print("local symbol == 0:", block_sym['local'] == {})
