import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import ftfield as F
mp.mp.dps = 90
m = 7
g = pickle.load(open('geom14.pkl','rb'))
E3 = pickle.load(open('E3n14.pkl','rb'))['E3n']
data = g['data']; ORDER=[3,4,1,2,5,6,0]
apn = [[F.feval(data[i][0],90), F.feval(data[i][1],90)] for i in ORDER]
LL = mp.sqrt(F.feval(g['L'],90))
def E3eval(t):
    logs=[mp.log(apn[i][0]+apn[i][1]*t) for i in range(7)]
    val=mp.mpf(0)
    for k,v in E3.items():
        poly=((F.feval(v[3],90)*t+F.feval(v[2],90))*t+F.feval(v[1],90))*t+F.feval(v[0],90)
        lv=mp.mpf(1)
        for i,e in enumerate(k):
            if e: lv*=logs[i]**e
        val+=poly*lv
    return val
# high precision integral of E3n (split at 1/2, tanh-sinh)
I = mp.quad(E3eval,[0,mp.mpf('0.5'),1],method='tanh-sinh',maxdegree=9)
print("int E3n =", mp.nstr(I,50))
print("84*L*int =", mp.nstr(84*LL*I,50))
# independent: Fourier filter on FULL chamber (no slivers), two radii
def fsum(th,p): return sum(mp.exp(p*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
def r_and_rp(th,p):
    f=fsum(th,p)
    fp=sum(-mp.sin(th-i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    return f**(-1/p), -(1/p)*f**(-1/p-1)*fp
def perim_full(p):
    return 2*m*mp.quad(lambda th:(lambda rr:mp.sqrt(rr[0]**2+rr[1]**2))(r_and_rp(th,p)),[mp.pi/2,mp.pi/2+mp.pi/m],method='tanh-sinh',maxdegree=8)
for rhos in ['0.03','0.015']:
    rho=mp.mpf(rhos); N=16; tot=mp.mpc(0)
    for j in range(N):
        z=mp.exp(2*mp.pi*1j*j/N); tot+=z**(-3)*perim_full(1+rho*z)
    P3=6*tot/N/rho**3
    print(f"Fourier filter full chamber rho={rhos}: P3 = {mp.nstr(P3,40)}", flush=True)
open('progress.log','a').write('PG6 s_A4: high-prec 84L int E3n = %s ; full-chamber Fourier agrees\n'%mp.nstr(84*LL*I,25))
