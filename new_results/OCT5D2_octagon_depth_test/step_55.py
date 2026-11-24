import os, sys, pickle, math, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import symengine as S
from symengine import Q, ZERO, ONE, P_LET, Q_LET, add, scale, shuffle, const_pow, obstruction, li_symbol, gsym
from fractions import Fraction as F

R2 = Q(0,1); B = Q(-1,1)
ZLIST = [('b',B),('-b',-B),('b^2',B*B),('-b^2',-B*B),('1/2',Q(F(1,2))),('1/sqrt2',R2/2),
         ('sqrt2',R2),('-sqrt2',-R2),('2-sqrt2',2-R2),('p',P_LET),('q',Q_LET),('-1',-ONE),('2',Q(2))]
t0=time.time()
cols = []  # (label, kind, symbol dict)
for name,z in ZLIST:
    cols.append((f"Li5({name})", 'Li5', li_symbol(5, z)))
for name,z in ZLIST:
    s4 = li_symbol(4, z)
    cols.append((f"Li4({name})*X", 'p41', shuffle(s4, const_pow(0,1))))
    cols.append((f"Li4({name})*Y", 'p41', shuffle(s4, const_pow(1,1))))
for name,z in ZLIST:
    s3 = li_symbol(3, z)
    cols.append((f"Li3({name})*X^2/2", 'p32', shuffle(s3, const_pow(0,2))))
    cols.append((f"Li3({name})*X*Y",   'p32', shuffle(s3, shuffle(const_pow(0,1),const_pow(1,1)))))
    cols.append((f"Li3({name})*Y^2/2", 'p32', shuffle(s3, const_pow(1,2))))
for name,z in ZLIST:
    s2 = li_symbol(2, z)
    cols.append((f"Li2({name})*X^3/6", 'p23', shuffle(s2, const_pow(0,3))))
    cols.append((f"Li2({name})*X^2Y/2", 'p23', shuffle(s2, shuffle(const_pow(0,2),const_pow(1,1)))))
    cols.append((f"Li2({name})*XY^2/2", 'p23', shuffle(s2, shuffle(const_pow(0,1),const_pow(1,2)))))
    cols.append((f"Li2({name})*Y^3/6", 'p23', shuffle(s2, const_pow(1,3))))
# log^5 terms
for i in range(6):
    cols.append((f"X^{i}Y^{5-i}/({i}!{5-i}!)", 'log5', shuffle(const_pow(0,i), const_pow(1,5-i))))
# Li32 atoms
def w32(x, y): return (ZERO, ZERO, ONE/x, ZERO, ONE/(x*y))
LI32 = [('(b,1)',(B,ONE)),('(1/sqrt2,1)',(R2/2,ONE)),('(-b,1)',(-B,ONE)),
        ('(b,-1)',(B,-ONE)),('(-b,-1)',(-B,-ONE)),('(1/sqrt2,-1)',(R2/2,-ONE)),
        ('(1-1/sqrt2,1)',(ONE-R2/2,ONE)),('(2-sqrt2,1)',(2-R2,ONE)),('(-b^2,1)',(-B*B,ONE))]
for name,(x,y) in LI32:
    cols.append((f"Li32{name}", 'li32', gsym(w32(x,y))))
print("columns:", len(cols), "built in %.1fs" % (time.time()-t0))

# matrix over words
raw = pickle.load(open("block_symbols.pkl","rb"))
def parse(s): return {tuple(map(int,w)): Q(F(p), F(q)) for w,(p,q) in s.items()}
W5sym = parse(raw['total'])
words = sorted(set(W5sym) | set().union(*[set(c[2]) for c in cols]))
print("symbol space dimension (words present):", len(words))
M = [[c[2].get(w, ZERO) for c in cols] for w in words]
bvec = [W5sym.get(w, ZERO) for w in words]
pickle.dump({'words': words, 'M': M, 'b': bvec, 'cols': [(l,k) for l,k,_ in cols]},
            open("symbol_fit.pkl","wb"))
print("matrix saved")
