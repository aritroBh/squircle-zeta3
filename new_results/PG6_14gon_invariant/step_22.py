import mpmath as mp
mp.mp.dps = 70
def fsum(th,p,m): return sum(mp.exp(p*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
def r_and_rp(th,p,m):
    f=fsum(th,p,m)
    fp=sum(-mp.sin(th-i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    return f**(-1/p), -(1/p)*f**(-1/p-1)*fp
def integrand(th,p,m):
    r,rp=r_and_rp(th,p,m); return mp.sqrt(r*r+rp*rp)
def perim(p,m,deg=10):
    hw=mp.pi/(2*m)
    return 2*m*mp.quad(lambda th: integrand(th,p,m),[-hw,0,hw],method='tanh-sinh',maxdegree=deg)
m=6
P1=perim(1,m)
print("P(1) 12gon =", mp.nstr(P1,30), " diff:", mp.nstr(P1-12*mp.sqrt(26-15*mp.sqrt(3)),3))
# Fourier filter for P3 (stays near p=1)
for rhos in ['0.03','0.015']:
    rho=mp.mpf(rhos); N=16; tot=mp.mpc(0)
    for j in range(N):
        z=mp.exp(2*mp.pi*1j*j/N); tot+=z**(-3)*perim(1+rho*z,m)
    P3=6*tot/N/rho**3
    print(f"12gon Fourier rho={rhos}: P3={mp.nstr(P3,35)}",flush=True)
known=mp.mpf('-5.20324042951766621029582970427243461666280533')
print("known =", mp.nstr(known,35))
# finite differences small h with split chamber
c3={-4:-mp.mpf(7)/240,-3:mp.mpf(3)/10,-2:-mp.mpf(169)/120,-1:mp.mpf(61)/30,1:-mp.mpf(61)/30,2:mp.mpf(169)/120,3:-mp.mpf(3)/10,4:mp.mpf(7)/240}
res={}
for hs in ['0.01','0.005','0.0025']:
    h=mp.mpf(hs); vals={k:perim(1+k*h,m) for k in range(-4,5)}
    res[hs]=sum(c3[k]*vals[k] for k in c3)/h**3
    print(f"12gon FD h={hs}: P3={mp.nstr(res[hs],35)}",flush=True)
print("FD vs known diffs:", [mp.nstr(abs(res[x]-known),3) for x in ['0.01','0.005','0.0025']])
open('progress.log','a').write('PG6 s_A6: 12-gon control (split chamber) FD+Fourier reproduce known P3 -5.2032404295176662\n')
