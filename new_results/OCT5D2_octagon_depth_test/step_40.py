import os, sys, time, pickle, math
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import symengine as S
from symengine import add, scale, shuffle, gsym, const_pow, obstruction
from fractions import Fraction as F
from functools import cache

@cache
def dual_full(n,a,b,c,d):
    out = {}
    for cp in range(c+1):
        for dp in range(d+1):
            binom = Q(math.comb(c,cp)*math.comb(d,dp))
            eX = (c-cp)+(d-dp)
            scal = binom * Q((-1)**(c-cp)) / Q(2**eX)
            parts = []
            if a: parts.append((ONE,)*a)
            if b: parts.append((ZERO,)*b)
            if cp: parts.append((Q_LET,)*cp)
            if dp: parts.append((P_LET,)*dp)
            sw = E.shuffle_words(parts)
            fac = scal * Q(math.factorial(a)*math.factorial(b)*math.factorial(cp)*math.factorial(dp))
            for k in range(n+1):
                sgn = Q((-1)**k) * Q(math.comb(n,k))
                for w, mult in sw.items():
                    for v, cv in E.IG(k, w).items():
                        key = (v, eX)
                        out[key] = out.get(key, ZERO) + sgn*fac*mult*cv
    return {k2: v for k2, v in out.items() if v != ZERO}

symcache = {}
def Gsym(w):
    if w not in symcache: symcache[w] = gsym(w)
    return symcache[w]
def dual_top_symbol(n,a,b,c,d):
    k = a+b+c+d
    out = {}
    for (v, eX), gam in dual_full(n,a,b,c,d).items():
        if len(v) + eX != k: continue
        s = Gsym(v)
        if eX: s = shuffle(s, const_pow(0, eX))
        out = add(out, scale(s, gam))
    return out

SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(SRC+"parsed_mons.pkl","rb"))
def qs2(p,q): return Q(F(p.numerator,p.denominator), F(q.numerator,q.denominator))

block_sym_d = {'local': {}, 'single': {}, 'mixed': {}}
t0=time.time()
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = 'mixed' if (c>0 and d>0) else ('single' if (c>0 or d>0) else 'local')
    for n,(p,q) in cu.items():
        coef = qs2(p,q)
        s = dual_top_symbol(n,a,b,c,d)
        if e: s = shuffle(s, const_pow(0, e))
        if f: s = shuffle(s, const_pow(1, f))
        block_sym_d[cls] = add(block_sym_d[cls], scale(s, coef))
block_sym_d = {k: scale(v, Q(960)) for k,v in block_sym_d.items()}
tot_d = add(block_sym_d['local'], block_sym_d['single'], block_sym_d['mixed'])
print("dual assembly in %.1fs" % (time.time()-t0))
for k in ['local','single','mixed']:
    print(k, "dual obstruction:", obstruction(block_sym_d[k]))
print("total dual obstruction:", obstruction(tot_d))
