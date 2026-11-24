import pickle, mpmath as mp, time
t0=time.time()
mp.mp.dps=400
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
def imli2(z): return mp.im(mp.polylog(2,z))
p2n = 1-mp.exp(1j*pi/4); q1n = r2+1j
base=[('z3',z3),('Lchi8',Lchi8),('pi3',pi**3),('r2*pi3',r2*pi**3),
      ('pi*G',pi*G),('pi*Cl1',pi*Cl1),
      ('pi*imLi2ib',pi*imli2(1j*b)),('pi*imLi2ia',pi*imli2(1j*a)),
      ('pi*imLi2p2',pi*imli2(p2n)),('pi*imLi2q1',pi*imli2(q1n)),
      ('z2*lam',z2*lam),('z2*L',z2*L),('z2*l3',z2*l3),
      ('r2*z3',r2*z3),('r2*Lchi8',r2*Lchi8),('r2*pi*Cl1',r2*pi*Cl1),('r2*pi*G',r2*pi*G)]
names=[nm for nm,_ in base]; vals=[v for _,v in base]
print("basis computed %.0fs"%(time.time()-t0))
K=pickle.load(open('K_V3p_400.pkl','rb'))
Kre=mp.mpf(K['re']); Kim=mp.mpf(K['im'])
mp.mp.dps=380
for part,nm in [(Kre,'Re K'),(Kim,'Im K')]:
    vec=[part]+vals
    rel=mp.pslq(vec, tol=mp.mpf('1e-340'), maxcoeff=10**20, maxsteps=2000)
    ok = rel and rel[0]!=0
    print(nm, "PSLQ:", "FOUND" if ok else "(none)")
    if ok:
        res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
        print("  residual:", mp.nstr(abs(res),6))
        print("  terms:", [(names[i-1] if i>0 else nm,r) for i,r in enumerate(rel) if r!=0])
print("time %.0fs"%(time.time()-t0))
