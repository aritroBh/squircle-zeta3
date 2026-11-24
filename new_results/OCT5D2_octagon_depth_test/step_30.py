import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, add, scale, shuffle, gsym, const_pow
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

cnt = 0
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a,b,c,d,e,f = mon
    if not (c==0 and d==0): continue
    for n,(p,q) in cu.items():
        sI = top_symbol(a,b,c,d,n)
        if sI:
            cnt += 1
            if cnt <= 3:
                print("nonzero local top_symbol:", mon, "n=",n)
                print("   sI =", sI)
print("local (mon,n) pairs with nonzero top symbol:", cnt)
