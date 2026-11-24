import sympy as sp, itertools
rt3=sp.sqrt(3)
s3=1/rt3
NB=[(1,0,0),(0,1,0),(0,0,1),(s3,s3,s3),(s3,s3,-s3),(s3,-s3,s3),(-s3,s3,s3)]
NB=[tuple(sp.nsimplify(x) for x in n) for n in NB]
a=(4*rt3-3)/13; c=(2*rt3-3)/2
verts=[];vlab=[]
for ax in range(3):
    for s in [1,-1]:
        v=[0,0,0];v[ax]=s*a;verts.append(tuple(v));vlab.append('A%d%s'%(ax,'+' if s>0 else '-'))
for zax in range(3):
    o=[i for i in range(3) if i!=zax]
    for s1 in [1,-1]:
        for s2 in [1,-1]:
            v=[0,0,0];v[o[0]]=s1*c;v[o[1]]=s2*c;verts.append(tuple(v));vlab.append('D%d%s%s'%(zax,'+' if s1>0 else '-','+' if s2>0 else '-'))
def avec(v): return [sp.simplify(sum(n[i]*v[i] for i in range(3))) for n in NB]
def dist2(x,y): return sp.simplify(sum((x[i]-y[i])**2 for i in range(3)))
L1sq=sp.Rational(2805,338)-sp.Rational(804,169)*rt3
L2sq=sp.Rational(21,2)-6*rt3
def show(i,j,tag):
    aa,ab=avec(verts[i]),avec(verts[j])
    print("%s %s->%s"%(tag,vlab[i],vlab[j]))
    print("   a_i=",[sp.nsimplify(x) for x in aa])
    print("   a_j=",[sp.nsimplify(x) for x in ab])
    for k in range(7):
        if aa[k]!=0 and ab[k]!=0:
            r=sp.simplify(ab[k]/aa[k])
            print("     term %d both: %s -> %s  ratio=%s = %s"%(k,sp.nsimplify(aa[k]),sp.nsimplify(ab[k]),sp.nsimplify(r),sp.radsimp(r)))
# one axis-diag edge
print("=== AXIS-DIAG edges (L1) ===")
cnt=0
for i in range(6):
    for j in range(6,18):
        if dist2(verts[i],verts[j])==L1sq:
            show(i,j,'L1'); cnt+=1
    if cnt>=2: break
print("\n=== DIAG-DIAG edges (L2) ===")
cnt=0
for i in range(6,18):
    for j in range(i+1,18):
        if dist2(verts[i],verts[j])==L2sq:
            show(i,j,'L2'); cnt+=1
    if cnt>=3: break
