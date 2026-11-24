from fractions import Fraction as Fr
# basis order: [1, ln2, ln2^2, ln2^3, pi^2, pi^2 ln2, pi^4, pi^4 ln2]
def Z(*args):
    v = [Fr(0)]*8
    for i,a in enumerate(args): v[i]=Fr(a)
    return tuple(v)
def zadd(x,y): return tuple(a+b for a,b in zip(x,y))
def zscale(s,x): return tuple(s*a for a in x)
def zmul(x,y):
    # multiply two z-values: ln2*ln2=ln2^2, ln2^2*ln2=ln2^3, pi^2*pi^2=pi^4, higher->0 (won't occur)
    # represent as (ln2pow, pi2pow) grid
    def split(z):
        d = {}
        idx = [(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(0,2),(1,2)]
        for i,a in enumerate(z): d[idx[i]]=a
        return d
    dx, dy = split(x), split(y)
    out = {}
    for (l1,p1),a in dx.items():
        for (l2,p2),b in dy.items():
            key=(l1+l2,p1+p2)
            out[key]=out.get(key,Fr(0))+a*b
    idx = [(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(0,2),(1,2)]
    return tuple(out.get(k,Fr(0)) for k in idx)

def bpv(n):
    if n==1: return Z(Fr(-1,2),Fr(1,2))  # (ln2-1)/2
    return Z(Fr((-1)**((n-1)//2), (n**2-1)*(2*((n-1)//4)+1)))
def bppv(n):
    if n==1: return Z(Fr(3,8),Fr(-1),Fr(1,4),0,Fr(5,96),0,Fr(-1,768),0)
    d = {3:(12,-16,3,128),5:(508,48,-45,3456),7:(-572,-96,63,13824),
         9:(2532,160,-225,64000),11:(-25592,-2160,2475,1296000),13:(1598344,75600,-143325,88905600)}
    a,b,c,dd = d[n]
    return Z(Fr(a,dd),Fr(b,dd),0,0,Fr(c,dd))
S13 = Z()
for n in [1,3,5,7,9,11,13]:
    S13 = zadd(S13, zmul(bpv(n), bppv(n)))
names = ["1","ln2","ln2^2","ln2^3","pi^2","pi^2 ln2","pi^4","pi^4 ln2"]
for nm,c in zip(names,S13):
    if c: print(nm, ":", c)
