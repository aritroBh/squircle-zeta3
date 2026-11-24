import os, sys, time, pickle
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
from fractions import Fraction as F
from linalg_tools import solve_linear
dd=pickle.load(open('det12_points.pkl','rb'))
orbits8=dd['orbits8']; mons8=dd['mons8']
data=pickle.load(open('eval12_0.pkl','rb'))+pickle.load(open('eval12_170.pkl','rb'))
N8=len(orbits8)
orbit_of8={}
for i,o in enumerate(orbits8):
    for e in o: orbit_of8[e]=i
def orbit_vals(a):
    vals=[0]*N8
    for e in mons8:
        pw=1
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        vals[orbit_of8[e]]+=pw
    return vals
def Lval(a):
    a1,a2,a3,a4,a5=a; return F(a2)*a3*a4*(a1+a2+a3+a4+a5)
def coords_of(d):
    g=d['g']; m0=d['m0']
    ft=[g[i]+m0[i] for i in range(12)]
    return {'A':ft[7],'B':ft[8],'CLASS':ft[9]+F(14,9)*ft[10],'Z':ft[11]}
from flint import nmod_mat
rows_all=[orbit_vals(d['pt']) for d in data]
p=(1<<61)-1
selidx=[]; cur=[]; r=0
for i,row in enumerate(rows_all):
    cur.append([int(x)%p for x in row])
    rr=nmod_mat(len(cur),N8,[x for rr2 in cur for x in rr2],p).rank()
    if rr>r:
        selidx.append(i); r=rr
        if r==N8: break
print("rank-255 subset from 340 pts:", len(selidx), "rank", r)
rows=[[F(x) for x in rows_all[i]] for i in selidx]
b=[coords_of(data[i])['A']/Lval(data[i]['pt']) for i in selidx]
t0=time.time()
sol=solve_linear(rows,b)
print("solve:", sol is not None, f"{time.time()-t0:.1f}s")
if sol:
    QA={}
    for oi,cf in zip(range(N8),sol):
        if cf:
            for e in orbits8[oi]:
                QA[e]=QA.get(e,F(0))+cf
    print("QA terms:", len(QA))
    def peval(Q,a):
        tot=F(0)
        for e,c in Q.items():
            pw=F(1)
            for x,ei in zip(a,e):
                if ei: pw*=x**ei
            tot+=c*pw
        return tot
    bad=0
    for d in data:
        D=coords_of(d)['A']-Lval(d['pt'])*peval(QA,d['pt'])
        if D!=0: bad+=1
    print("D_A==0 on all 340 pts:", bad==0, "bad=",bad)
    import sympy as sp
    print("QA(1,1,1,1,1) =", peval(QA,(1,1,1,1,1)))
