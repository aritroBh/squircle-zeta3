import pickle, mpmath as mp, time
from q1s_exact import femb
t0=time.time()
mp.mp.dps=400
sols=pickle.load(open('sunit_li3.pkl','rb'))
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
alpha=mp.atan(1/r2)
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; z2=pi**2/6
def re3(i): return mp.re(mp.polylog(3,femb(sols[i]['z'],1)))
def im3(i): return mp.im(mp.polylog(3,femb(sols[i]['z'],1)))
R=pickle.load(open('Rreal_V3p_400.pkl','rb'))
Rre=mp.mpf(R['Rre'])
mp.mp.dps=380
def try_pslq(basis, label, tol='1e-320', mc=10**24):
    names2=[nm for nm,_ in basis]; vals=[mp.mpf(v) for _,v in basis]
    vec=[Rre]+vals
    rel=mp.pslq(vec, tol=mp.mpf(tol), maxcoeff=mc, maxsteps=3000)
    ok = rel and rel[0]!=0
    print(label, ":", "FOUND" if ok else "(none)")
    if ok:
        res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
        print("  residual:", mp.nstr(abs(res),6))
        print("  terms:", [(names2[i-1] if i>0 else 'R',r) for i,r in enumerate(rel) if r!=0])
    return ok
core=[('z3',z3),('Lchi8',Lchi8),('pi3',pi**3),('r2*pi3',r2*pi**3),('pi*G',pi*G),
      ('pi*Cl1',mp.pi*mp.clsin(2,mp.pi/4)),
      ('z2*lam',z2*lam),('z2*L',z2*L),('z2*l3',z2*l3),('z2*alpha',z2*alpha)]
try_pslq(core, "core")
core2=core+[('Cl3(pi/4)',mp.re(mp.polylog(3,mp.exp(1j*pi/4)))),('Cl3(3pi/4)',mp.re(mp.polylog(3,mp.exp(3j*pi/4))))]
try_pslq(core2, "core+Cl3")
tors=[243,244,245,246,247,248,249,289,292,376,377,399,401,479]
core3=core2+[(f'reLi3(z{i})',re3(i)) for i in tors]+[(f'imLi3(z{i})',im3(i)) for i in tors]
core3=[(nm,v) for nm,v in core3 if abs(mp.mpf(v))>mp.mpf('1e-300')]
try_pslq(core3, "core+Cl3+torsLi3")
core4=core3+[('pi*Cl2(pi/8)',pi*mp.clsin(2,pi/8)),('pi*Cl2(3pi/8)',pi*mp.clsin(2,3*pi/8))]
try_pslq(core4, "core+Cl3+tors+Cl2(pi/8,3pi/8)")
print("time %.0fs"%(time.time()-t0))
