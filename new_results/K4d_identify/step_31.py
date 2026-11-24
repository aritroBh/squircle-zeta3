import sympy as sp
import itertools, collections
rt2 = sp.sqrt(2)
# exact vertices
val4 = []
for ax in range(3):
    for s in [1,-1]:
        v=[0,0,0]; v[ax]=s*rt2/4; val4.append(tuple(v))
val6 = []
for signs in itertools.product([1,-1],repeat=3):
    val6.append(tuple(s*rt2/6 for s in signs))
verts = val4+val6
labels = ['A%d'%i for i in range(6)]+['B%d'%i for i in range(8)]
# normals (6): permutations of (1,1,0)/sqrt2 with independent signs
normals=[]
for perm in [(0,1,2),(0,2,1),(1,2,0)]:
    pass
# build as in NA: (1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1) /sqrt2
NA=[(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
NA=[tuple(sp.Rational(x)/rt2 for x in n) for n in NA]
def aval(n,v): return sp.simplify(sum(n[i]*v[i] for i in range(3)))
def dist2(a,b): return sp.simplify(sum((a[i]-b[i])**2 for i in range(3)))
# edges by distance
d2vals = collections.Counter()
pairs=[]
for i,j in itertools.combinations(range(14),2):
    d2=dist2(verts[i],verts[j])
    d2vals[d2]+=1
print("pairwise dist^2 histogram:", {str(k):v for k,v in d2vals.items()})
# edges = two shortest distances
L_ab2 = sp.Rational(1,8)  # (1/(2sqrt2))^2
L_bb2 = sp.Rational(2,9)  # (sqrt2/3)^2
edges=[]
for i,j in itertools.combinations(range(14),2):
    d2=dist2(verts[i],verts[j])
    if d2 in (L_ab2, L_bb2):
        edges.append((i,j,d2))
print("num edges:", len(edges), " apex-base:", sum(1 for e in edges if e[2]==L_ab2), " base-base:", sum(1 for e in edges if e[2]==L_bb2))
