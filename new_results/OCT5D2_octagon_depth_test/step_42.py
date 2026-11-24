import os, sys, time, pickle, math
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import symengine as S
from symengine import add, scale, shuffle, gsym, const_pow, obstruction
from fractions import Fraction as F

SRC = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(SRC+"parsed_mons.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
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

def sub(mon):
    a,b,c,d,e,f = mon
    if c>0 and d>0: return 'mixed'
    if c>0: return 'single-p'
    if d>0: return 'single-q'
    return 'local'
sub_sym = {}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    cls = sub(mon)
    for n,(p,q) in cu.items():
        coef = qs2(p,q)
        s = top_symbol(a,b,c,d,n)
        if e: s = shuffle(s, const_pow(0, e))
        if f: s = shuffle(s, const_pow(1, f))
        sub_sym[cls] = add(sub_sym.get(cls,{}), scale(s, coef))
for k,v in sub_sym.items():
    print(f"{k:9s} obstruction:", obstruction(scale(v, Q(960))))
