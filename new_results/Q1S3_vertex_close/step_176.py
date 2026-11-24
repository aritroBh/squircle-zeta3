import pickle, mpmath as mp, time
from q1s_exact import femb
mp.mp.dps=300
sols=pickle.load(open('sunit_li3.pkl','rb'))
def re3(i): return mp.re(mp.polylog(3,femb(sols[i]['z'],1)))
def im3(i): return mp.im(mp.polylog(3,femb(sols[i]['z'],1)))
r2=mp.sqrt(2); pi=mp.pi
z3=mp.zeta(3)
Lchi8=(mp.zeta(3,mp.mpf(1)/8)-mp.zeta(3,mp.mpf(3)/8)-mp.zeta(3,mp.mpf(5)/8)+mp.zeta(3,mp.mpf(7)/8))/512
R=pickle.load(open('Rreal2_V3p_1250.pkl','rb'))
Rre=mp.mpf(R['Rre'])
tors=[243,244,245,246,247,248,249,289,292,376,377,399,401,479]
base=[('z3',z3),('Lchi8',Lchi8),('pi3',pi**3),('r2*pi3',r2*pi**3)]
for i in tors:
    base.append(('reLi3(z%d)'%i, re3(i)))
    base.append(('imLi3(z%d)'%i, im3(i)))
base=[(nm,mp.mpf(v)) for nm,v in base if abs(mp.mpf(v))>mp.mpf('1e-280')]
names2=[nm for nm,_ in base]; vals=[v for _,v in base]
print("basis:", len(vals))
# try mp.pslq with relaxed tol and more steps
vec=[Rre]+vals
for tol,mc,ms in [(mp.mpf('1e-150'),10**6,5000),(mp.mpf('1e-200'),10**3,8000)]:
    rel=mp.pslq(vec, tol=tol, maxcoeff=mc, maxsteps=ms)
    ok = rel and rel[0]!=0
    print("tol",mp.nstr(tol,3),"mc",mc,":", "FOUND" if ok else "(none)")
    if ok:
        res=sum(mp.mpf(r)*v for r,v in zip(rel,vec))
        print("  residual:", mp.nstr(abs(res),6))
        print("  terms:", [(names2[i-1] if i>0 else 'R',r) for i,r in enumerate(rel) if r!=0])
        break
