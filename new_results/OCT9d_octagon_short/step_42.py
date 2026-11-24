import sys, time, itertools, pickle, random
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()
Gsym = {g: col_Li3(sp.Rational(g), SS235) for g in Gargs}
sel = ['1/3', '2/3', '2/5', '-1/5', '1/10', '-5/3', '5', '3/5']
selq = [sp.Rational(g) for g in sel]
G8 = sp.Matrix.hstack(*[Gsym[g] for g in selq])

# candidate product columns: nice order
EP = set()
for e in (1,2,3):
    for z in edge_points(e): EP.add(z)
ORB = set()
for z in list(EP): ORB |= orbit(z)
ORBSM = set()
for z in list(SMALL235): ORBSM |= orbit(z)
def nice(g):
    g = Fr(g)
    if g in EP: return 0
    if g in ORB: return 1
    if g in SMALL235: return 2
    if g in ORBSM: return 3
    return 4
wpairs = [(w,p) for w in sorted(Gargs, key=lambda g:(nice(g), max(abs(Fr(g).numerator),abs(Fr(g).denominator)), str(g))) for p in (2,3,5)]
pcolmap = {}
def pcol(w,p):
    if (w,p) not in pcolmap: pcolmap[(w,p)] = col_prod(sp.Rational(w), p, SS235)
    return pcolmap[(w,p)]

cur = []
curmat = sp.zeros(27,8)  # G8 included implicitly
curcols = sp.Matrix.hstack(G8, sp.zeros(27,0))
for (w,p) in wpairs:
    c = pcol(w,p)
    if sp.Matrix.hstack(curcols, c).rank() > curcols.shape[1]:
        cur.append((w,p)); curcols = sp.Matrix.hstack(curcols, c)
    if len(cur) == 19: break
def fulldet(sel19):
    M = sp.Matrix.hstack(G8, *[pcol(w,p) for w,p in sel19])
    return abs(M.det())
d = fulldet(cur)
print("start det:", d, "=", sp.factorint(d), f"({time.time()-t0:.1f}s)")
improved = True
it=0
while improved and it < 12:
    improved = False; it += 1
    for i in range(19):
        for (w,p) in wpairs:
            if (w,p) in cur: continue
            sel2 = cur[:i]+[(w,p)]+cur[i+1:]
            d2 = fulldet(sel2)
            if d2 != 0 and d2 < d:
                cur, d = sel2, d2; improved = True
        print(f"  pass {it} i={i}: det {d} = {sp.factorint(d)} ({time.time()-t0:.1f}s)")
print("final PB:", [(str(w),p) for w,p in cur])
print("final det:", d, sp.factorint(d))
pickle.dump({'PB19': [(str(w),p) for w,p in cur], 'det': str(d)}, open("PB19.pkl","wb"))
