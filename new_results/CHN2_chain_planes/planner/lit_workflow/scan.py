import sys, pickle; sys.dont_write_bytecode=True
import numpy as np, sympy as sp
a1,a2,a3,a4=sp.symbols('a1 a2 a3 a4')
K53,K73,K75=pickle.load(open('K.pkl','rb'))
f53=sp.lambdify((a1,a2,a3,a4),K53,'numpy'); f73=sp.lambdify((a1,a2,a3,a4),K73,'numpy'); f75=sp.lambdify((a1,a2,a3,a4),K75,'numpy')
P=sp.Poly(K53.subs({a4:1}),a1)
coeffs=[sp.lambdify((a2,a3),c,'numpy') for c in P.all_coeffs()]
g=np.geomspace(0.05,20,120)
pts=[]
for x2 in g:
    for x3 in g:
        cs=[float(c(x2,x3)) for c in coeffs]
        for r in np.roots(cs):
            if abs(r.imag)<1e-9 and r.real>0:
                pts.append((r.real,x2,x3,1.0))
pts=np.array(pts)
print('positive points on {K53=0} (a4=1 chart):',len(pts))
v73=np.array([f73(*p) for p in pts]); v75=np.array([f75(*p) for p in pts])
# normalise by |a|^5, |a|^7
nrm=np.array([np.linalg.norm(p) for p in pts])
s73=v73/nrm**5; s75=v75/nrm**7
# distance to Pi1 (a1=a2,a3=a4) and Pi2 (a1=a3,a2=a4) in the a4=1 chart
d1=np.sqrt((pts[:,0]-pts[:,1])**2+(pts[:,2]-1)**2); d2=np.sqrt((pts[:,0]-pts[:,2])**2+(pts[:,1]-1)**2)
dist=np.minimum(d1,d2)
far=dist>0.05
print('points farther than 0.05 from Pi1,Pi2:',far.sum())
print('on those: min |K73|/|a|^5 = %.3e ; min |K75|/|a|^7 = %.3e'%(np.min(np.abs(s73[far])),np.min(np.abs(s75[far]))))
# sign structure of K73 on the surface away from loci
print('K73 sign on surface (away from loci): +%d / -%d'%((s73[far]>0).sum(),(s73[far]<0).sum()))
idx=np.argsort(np.abs(s73[far]))[:5]
print('closest-to-zero K73 points away from loci (a1,a2,a3,a4) and dist:')
for i in idx: print('   ',np.round(pts[far][i],4),'K73n=%.3e'%s73[far][i],'K75n=%.3e'%s75[far][i],'dist=%.3f'%dist[far][i])
print('--- stricter: dist>0.3 from both lines')
far=dist>0.3
print('count',far.sum(),'min |K73n| %.3e  max |K73n| %.3e'%(np.min(np.abs(s73[far])),np.max(np.abs(s73[far]))))
print('sign +%d / -%d'%((s73[far]>0).sum(),(s73[far]<0).sum()))
# sign vs which side of the loci: a1<a3 vs a1>a3 etc
side=np.sign((pts[:,0]-pts[:,2]))*np.sign((pts[:,0]-pts[:,1]))
for sgn in (1,-1):
    m=far&(side==sgn)
    print('side',sgn,'count',m.sum(),'K73 sign + %d - %d'%((s73[m]>0).sum(),(s73[m]<0).sum()))
