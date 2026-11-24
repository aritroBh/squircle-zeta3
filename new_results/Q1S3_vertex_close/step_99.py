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
r2=mp.sqrt(2); a=1+r2; L=mp.log(a); lam=mp.log(2); l3=mp.log(3); pi=mp.pi
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
Lchim8=(mp.zeta(3,mp.mpf(1)/8)+mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)-mp.zeta(3,mp.mpf(7)/8))/512
beta3=(mp.zeta(3,mp.mpf(1)/4)-mp.zeta(3,mp.mpf(3)/4))/64
G=mp.catalan; Cl1=mp.clsin(2,mp.pi/4); z2=pi**2/6
base=[('z3',z3),('Lchi8',Lchi8),('Lchim8',Lchim8),('beta3',beta3),('pi3',pi**3),
      ('pi*G',pi*G),('pi*Cl1',pi*Cl1),
      ('z2*lam',z2*lam),('z2*L',z2*L),('z2*l3',z2*l3)]
names=[nm for nm,_ in base]; vals=[v for _,v in base]
for part,nm in [(mp.re(K),'Re K'),(mp.im(K),'Im K')]:
    vec=[part]+vals
    rel=mp.pslq(vec, tol=mp.mpf('1e-88'), maxcoeff=10**12, maxsteps=1000)
    print(nm, "PSLQ:", rel)
    if rel:
        res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
        print("  residual:", mp.nstr(abs(res),6), " terms:", [(names[i-1] if i>0 else nm,r) for i,r in enumerate(rel) if r!=0])
