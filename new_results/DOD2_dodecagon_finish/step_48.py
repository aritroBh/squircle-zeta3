import sympy as sp
r3 = sp.sqrt(3)
Xs = {1: 1-2*r3/3, 2: 1-r3, 5: 1-r3/3, 6: 1-r3/2}
IDX=[1,2,5,6]
gamS = {}
for i in IDX:
    for j in IDX:
        if i==j: continue
        gamS[(i,j)] = sp.simplify(Xs[j]/(Xs[j]-Xs[i]))
triples = []
for i in IDX:
    for j in IDX:
        for k in IDX:
            if len({i,j,k})==3 and j<k:
                triples.append((i,j,k))
gdset = set(gamS.values())
# add deltas
for (i,j) in list(gamS.keys()):
    gdset.add(sp.simplify(gamS[(i,j)]*(1-Xs[i])))
eset = set()
for g_ in list(gdset):
    eset.add(sp.simplify((g_-1)/g_))
xset = set(Xs.values()) | set(sp.simplify(1-Xs[i]) for i in IDX) | {r3/2, r3/3}
cr = {}
for (i,j,k) in triples:
    r_ = sp.simplify(gamS[(i,k)]/gamS[(i,j)])
    cr[(i,j,k)] = r_
print("cross-ratio points per triple:")
allcr = set()
for (i,j,k) in triples:
    r_ = cr[(i,j,k)]
    allcr.add(r_)
    tag = []
    if r_ in gdset: tag.append("gd")
    if r_ in eset: tag.append("e")
    if r_ in xset: tag.append("x")
    print((i,j,k), "cr =", r_, "~", sp.N(r_,6), "sets:", tag or "NEW")
print("distinct cross-ratios:", len(allcr))
for r_ in sorted(allcr, key=lambda e: float(sp.N(e,10))):
    tag=[]
    if r_ in gdset: tag.append("gd")
    if r_ in eset: tag.append("e")
    if r_ in xset: tag.append("x")
    print("   ", r_, "~", sp.N(r_,6), tag or "NEW")
