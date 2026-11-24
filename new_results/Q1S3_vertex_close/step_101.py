import pickle, mpmath as mp, json
from q1s_exact import femb, BASISEL
mp.mp.dps=100
sol = pickle.load(open('q3_solve.pkl','rb'))
sols = pickle.load(open('sunit_li3.pkl','rb'))
def ReLi3(i): return mp.re(mp.polylog(3, femb(sols[i]['z'],1)))
def logu(u): return mp.log(femb(BASISEL[u],1))
def Gvw(v,w):
    vn=femb(BASISEL[v],1); wn=femb(BASISEL[w],1)
    return mp.quad(lambda t: mp.log(1-t/wn)/(t-vn), [0, mp.mpf('0.5'), 1])
li3part = sum(mp.mpf(float(c))*2*ReLi3(i) for (i,j),c in sol['V3p']['li3'])
decpart = sum(mp.mpf(float(c))*logu(u)*Gvw(v,w) for (u,v,w),c in sol['V3p']['dec'])
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
K = ref3 - li3part - decpart
r2=mp.sqrt(2); a=1+r2; b=r2-1; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
Lchim8=(mp.zeta(3,mp.mpf(1)/8)+mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)-mp.zeta(3,mp.mpf(7)/8))/512
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
def imli2(z): return mp.im(mp.polylog(2,z))
p2n = 1-mp.exp(1j*pi/4); q1n = r2+1j; q2n = r2-1j
base=[('z3',z3),('Lchi8',Lchi8),('Lchim8',Lchim8),('pi3',pi**3),
      ('pi*G',pi*G),('pi*Cl1',pi*Cl1),
      ('z2*lam',z2*lam),('z2*L',z2*L),('z2*l3',z2*l3),
      ('pi*imLi2ib',pi*imli2(1j*b)),('pi*imLi2ia',pi*imli2(1j*a)),
      ('pi*imLi2p2',pi*imli2(p2n)),('pi*imLi2q1',pi*imli2(q1n)),
      ('pi*lam2',pi*lam**2),('pi*L2',pi*L**2),('pi*lam*L',pi*lam*L),
      ('z2',z2)]
names=[nm for nm,_ in base]; vals=[v for _,v in base]
for part,nm in [(mp.re(K),'Re K'),(mp.im(K),'Im K')]:
    vec=[part]+vals
    rel=mp.pslq(vec, tol=mp.mpf('1e-80'), maxcoeff=10**12, maxsteps=2000)
    ok = rel and rel[0]!=0
    print(nm, "PSLQ:", rel if ok else "(none)")
    if ok:
        res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
        print("  residual:", mp.nstr(abs(res),6), " terms:", [(names[i-1] if i>0 else nm,r) for i,r in enumerate(rel) if r!=0])
