import sys, mpmath as mp, pickle
sys.path.insert(0,'.')
mp.mp.dps=60
import engine14 as E
# numerical values of atoms
Xv={i: E.fev(E.X[i],60) for i in E.IDX}
l2v=mp.log(2); lav=mp.log(2*mp.cos(2*mp.pi/7)); luv=mp.log(2*mp.cos(mp.pi/7))
z2v=mp.pi**2/6; z3v=mp.zeta(3)
w={i:(lambda tt,i=i: mp.log(1-Xv[i]*tt)) for i in E.IDX}
u={i:(lambda tt,i=i: 1-Xv[i]*tt) for i in E.IDX}
def Li21(x,y):  # Goncharov prefix f(x,y)=sum_{0<k<m} x^m y^k/(m^2 k)
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1)
    for m in range(1,6000):
        xp*=x
        if m>=2: tot+=xp*ps/m**2
        yp*=y; ps+=yp/m
    return tot
def Aser(x):
    tot=mp.mpf(0); xp=mp.mpf(1); H=mp.mpf(0)
    for n in range(1,8000): xp*=x; H+=mp.mpf(1)/n; tot+=H*xp/n**2
    return tot
NUM={'l2':l2v,'la':lav,'lu':luv,'z2':z2v,'z3':z3v}
for i in E.IDX:
    NUM['L2x%d'%i]=mp.polylog(2,Xv[i]); NUM['L3x%d'%i]=mp.polylog(3,Xv[i]); NUM['Ax%d'%i]=Aser(Xv[i])
for (i,j) in E.Lam: NUM['Lam%d_%d'%(i,j)]=Li21(Xv[i],Xv[j]/Xv[i])
for (i,j) in E.T2s:  NUM['T2_%d_%d'%(i,j)]=mp.quad(lambda tt: w[j](tt)/u[i](tt),[0,1])
for (i,j) in E.T3ps: NUM['T3p_%d_%d'%(i,j)]=mp.quad(lambda tt: w[j](tt)**2/u[i](tt),[0,1])
for (i,j,k) in E.T3a: NUM['T3_%d_%d%d'%(i,j,k)]=mp.quad(lambda tt: w[j](tt)*w[k](tt)/u[i](tt),[0,1])
pickle.dump({k:str(v) for k,v in NUM.items()}, open('atoms_num.pkl','wb'))
def e_eval(P):
    tot=mp.mpf(0)
    for mon,c in P.items():
        v=E.fev(c,60)
        for (s,ex) in mon:
            v*=NUM[s]**ex
        tot+=v
    return tot
# test integral types vs quad
def quad_int(n,a,b,wc):
    f=lambda tt: tt**n*(mp.log(tt))**a*(mp.log(1-tt))**b*mp.prod([w[idx](tt)**c for idx,c in zip(E.IDX,wc) if c>0])
    return mp.quad(f,[0,1])
tests=[(0,0,0,(0,0,0,0)),(2,1,1,(0,0,0,0)),(0,0,1,(0,0,0,0)),
       (0,0,0,(1,0,0,0)),(1,1,0,(2,0,0,0)),(2,0,1,(0,0,1,0)),
       (0,0,0,(1,1,0,0)),(1,1,0,(1,1,0,0)),(0,0,0,(2,1,0,0)),(2,0,0,(1,2,0,0)),
       (0,0,0,(1,1,1,0)),(1,0,0,(1,1,1,0)),(3,0,0,(0,1,1,1))]
worst=mp.mpf(0)
for typ in tests:
    n,a,b,wc=typ
    got=e_eval(E.INT(n,a,b,wc)); tru=quad_int(n,a,b,wc)
    d=abs(got-tru); worst=max(worst,d)
    print(typ, "got", mp.nstr(got,15), "tru", mp.nstr(tru,15), "diff", mp.nstr(d,2))
print("WORST:", mp.nstr(worst,3))
open('progress.log','a').write('PG6 s_B5: engine14 validated vs quad, worst diff %s over %d types\n'%(mp.nstr(worst,2),len(tests)))
