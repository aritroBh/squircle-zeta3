import mpmath as mp
mp.mp.dps = 60
def fsum(th,p,m): return sum(mp.exp(p*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
def r_and_rp(th,p,m):
    f=fsum(th,p,m)
    fp=sum(-mp.sin(th-i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    return f**(-1/p), -(1/p)*f**(-1/p-1)*fp
def integrand(th,p,m):
    r,rp=r_and_rp(th,p,m); return mp.sqrt(r*r+rp*rp)
m=6
# check P(1) vs 12 sqrt(26-15 sqrt3)
P1 = 2*m*mp.quad(lambda th: integrand(th,1,m),[-mp.pi/(2*m),mp.pi/(2*m)],method='tanh-sinh',maxdegree=8)
print("P(1) 12gon =", mp.nstr(P1,30))
print("12 sqrt(26-15 sqrt3) =", mp.nstr(12*mp.sqrt(26-15*mp.sqrt(3)),30))
# check perim at the stencil points
for k in range(-4,5):
    for hs in ['0.02','0.0025']:
        h=mp.mpf(hs); p=1+k*h
        v=2*m*mp.quad(lambda th: integrand(th,p,m),[-mp.pi/(2*m),mp.pi/(2*m)],method='tanh-sinh',maxdegree=8)
        if abs(v)>3: print(f"  k={k} h={hs} p={mp.nstr(p,6)}: perim={mp.nstr(v,12)}  <-- LARGE")
print("done scan (only large values printed)")
# sample integrand across the chamber at p=1
for ths in ['-0.26','-0.13','0','0.13','0.26']:
    th=mp.mpf(ths)
    print("integrand(th=%s,p=1)=%s"%(ths, mp.nstr(integrand(th,1,m),15)))
print("pi/(2m)=", mp.nstr(mp.pi/12,10))
