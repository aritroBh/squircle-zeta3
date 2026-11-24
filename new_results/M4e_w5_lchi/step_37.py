import json, math, pickle, time
from fractions import Fraction as F
from math import gcd, lcm
import flint
from flint import arb, ctx
import mpmath as mp
import sys
M4E = '<project>/research/08_theory/M4e_w5_lchi'
sys.path.insert(0, M4E)
import w5_symbol as W
from w5_search import evaluate, mpl, numeric
ZERO, ONE, Q, R, A, B = W.ZERO, W.ONE, W.Q, W.R, W.A, W.B

digits = 500   # diagnostic
ctx.prec = math.ceil((digits + 140) * math.log2(10))
mp.mp.dps = digits + 120

def qnum(c):
    return mp.mpf(c.p.numerator)/c.p.denominator + mp.sqrt(2)*mp.mpf(c.q.numerator)/c.q.denominator
def prim(vec):
    D = lcm(*[c.denominator for c in vec]); ints=[int(c*D) for c in vec]; g=0
    for i in ints: g=gcd(g,abs(i))
    return [i//g for i in ints]

sym = json.loads(open(M4E+'/nice_symbol.json').read())
names = sym['names']; xp=[F(s) for s in sym['xp']]; xq=[F(s) for s in sym['xq']]
K=[prim([F(c) for c in v]) for v in sym['kernel']]
red = flint.fmpz_mat(K).lll(gram='exact')
K=[[int(red[i,j]) for j in range(red.ncols())] for i in range(red.nrows())]
recipes = pickle.loads(open(M4E+'/nice_atoms.pkl','rb').read())
values=[evaluate(r, digits+60) for _,r in recipes]
vals=[mp.mpf(v.mid().str(digits+110, radius=False)) for v in values]
r2=mp.sqrt(2)
cert=json.loads(open(M4E+'/cobracket.json').read())
fixedv=mp.mpf(0)
for term in cert['correction']:
    x,y,c=(Q(*map(F,term[k])) for k in ['x','y','coefficient'])
    v=evaluate(('mpl',3,2,x,y),digits+60)
    fixedv+=qnum(c)*mp.mpf(v.mid().str(digits+110, radius=False))
saved=json.loads(open(M4E+'/W5_5400.json').read())
pval=mp.fsum((mp.mpf(x.numerator)/x.denominator+r2*mp.mpf(q.numerator)/q.denominator)*v for x,q,v in zip(xp,xq,vals))
target=mp.mpf(saved['value'])-fixedv-pval
kv=[mp.fsum(F(c)*v for c,v in zip(krow,vals)) for krow in K]
Lchi=mp.mpf(json.loads(open(M4E+'/Lchi8_5.json').read())['value'])
colvals=[target]+kv+[r2*v for v in kv]+[Lchi, r2*Lchi]
# drop the known dependent columns kv16,kv28,kv29 & s2 partners (indices 17,29,30,52,64,65)
drop={17,29,30,52,64,65}
colvals=[v for j,v in enumerate(colvals) if j not in drop]
n=len(colvals)
C=mp.mpf(10)**digits
rows=[[int(i==j) for j in range(n)]+[int(mp.nint(C*v))] for i,v in enumerate(colvals)]
t0=time.monotonic()
redm=flint.fmpz_mat(rows).lll(gram='exact')
print(f"LLL {n}x{n} at {digits}: {time.monotonic()-t0:.1f}s")
out=[]
for i in range(n):
    vec=[int(redm[i,j]) for j in range(n)]
    residual=mp.fsum(c*v for c,v in zip(vec,colvals))
    out.append((vec[0], max(map(abs,vec)), residual, vec[-2], vec[-1]))
out.sort(key=lambda t: abs(t[2]) if t[2]!=0 else mp.inf)
for r in out[:12]:
    print(f"vec0={r[0]} maxcoef={r[1]} resid={mp.nstr(r[2],4)} Lcoefs=({r[3]},{r[4]})")
