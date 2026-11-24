import sympy as sp, itertools, collections
rt3=sp.sqrt(3)
# Model B normals
s3=1/rt3
NB=[(1,0,0),(0,1,0),(0,0,1),(s3,s3,s3),(s3,s3,-s3),(s3,-s3,s3),(-s3,s3,s3)]
NB=[tuple(sp.nsimplify(x) for x in n) for n in NB]
# Vertices: 6 axis at +-a e_i, a=(4rt3-3)/13 ; 12 at (0,+-c,+-c) perms, c=(2rt3-3)/2
a=(4*rt3-3)/13
c=(2*rt3-3)/2
verts=[]
vlab=[]
for ax in range(3):
    for s in [1,-1]:
        v=[0,0,0]; v[ax]=s*a; verts.append(tuple(v)); vlab.append('ax%d%s'%(ax,'+' if s>0 else '-'))
for zero_ax in range(3):
    others=[i for i in range(3) if i!=zero_ax]
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            v=[0,0,0]; v[others[0]]=s1*c; v[others[1]]=s2*c; verts.append(tuple(v)); vlab.append('d%d%s%s'%(zero_ax,'+' if s1>0 else '-','+' if s2>0 else '-'))
print("num vertices:", len(verts))
def avec(v): return [sp.simplify(sum(n[i]*v[i] for i in range(3))) for n in NB]
# check sum|a|=1 and print a couple
for i in [0,6]:
    av=avec(verts[i])
    print(vlab[i], verts[i], "a=", [sp.nsimplify(x) for x in av], "sum|.|=",sp.simplify(sum(abs(x) for x in av)))
# edges by distance
def dist2(a,b): return sp.simplify(sum((a[i]-b[i])**2 for i in range(3)))
d2c=collections.Counter()
for i,j in itertools.combinations(range(len(verts)),2):
    d2c[sp.nsimplify(dist2(verts[i],verts[j]))]+=1
print("dist^2 hist:", {str(k):v for k,v in sorted(d2c.items(), key=lambda kv: float(kv[0]))})
