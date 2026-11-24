import sys, mpmath as mp, pickle
sys.path.insert(0,'.')
mp.mp.dps=60
import engine14 as E
Xv={i: E.fev(E.X[i],60) for i in E.IDX}
w={i:(lambda tt,i=i: mp.log(1-Xv[i]*tt)) for i in E.IDX}
NUM={k:mp.mpf(v) for k,v in pickle.load(open('atoms_num.pkl','rb')).items()}
def e_eval(P):
    tot=mp.mpf(0)
    for mon,c in P.items():
        v=E.fev(c,60)
        for (s,ex) in mon: v*=NUM[s]**ex
        tot+=v
    return tot
def quad_int(n,a,b,wc):
    def f(tt):
        r=tt**n*(mp.log(tt))**a*(mp.log(1-tt))**b
        for idx,c in zip(E.IDX,wc):
            if c>0: r*= w[idx](tt)**c
        return r
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
    print(typ, "diff", mp.nstr(d,2))
print("WORST:", mp.nstr(worst,3))
open('progress.log','a').write('PG6 s_B5: engine14 validated vs quad, worst diff %s over %d types\n'%(mp.nstr(worst,2),len(tests)))
