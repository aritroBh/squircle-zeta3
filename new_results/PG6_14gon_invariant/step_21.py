import mpmath as mp
mp.mp.dps = 60
def fsum(th,p,m): return sum(mp.exp(p*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
def r_and_rp(th,p,m):
    f=fsum(th,p,m)
    fp=sum(-mp.sin(th-i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    return f**(-1/p), -(1/p)*f**(-1/p-1)*fp
def integrand(th,p,m):
    r,rp=r_and_rp(th,p,m); return mp.sqrt(r*r+rp*rp)
def perim(p,m,deg):
    return 2*m*mp.quad(lambda th: integrand(th,p,m),[-mp.pi/(2*m),mp.pi/(2*m)],method='tanh-sinh',maxdegree=deg)
m=6
for deg in [8,10,12]:
    P1=perim(1,m,deg)
    print(f"deg={deg}: P(1)={mp.nstr(P1,25)} diff={mp.nstr(P1-12*mp.sqrt(26-15*mp.sqrt(3)),3)}")
# sample integrand at p=0.92 to find the sharp feature
print("--- integrand profile at p=0.92 ---")
for i in range(-10,11):
    th = i/10*mp.pi/12
    print(f"  th={mp.nstr(th,5)}: {mp.nstr(integrand(th,mp.mpf('0.92'),m),10)}")
