import mpmath as mp
mp.mp.dps = 80
def fsum(th,p,m): return sum(mp.exp(p*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
def r_and_rp(th,p,m):
    f=fsum(th,p,m)
    fp=sum(-mp.sin(th-i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    return f**(-1/p), -(1/p)*f**(-1/p-1)*fp
def integrand(th,p,m):
    r,rp=r_and_rp(th,p,m); return mp.sqrt(r*r+rp*rp)
def perim(p,m):
    return 2*m*mp.quad(lambda th: integrand(th,p,m),[-mp.pi/(2*m),mp.pi/(2*m)],method='tanh-sinh',maxdegree=8)
m=6
c3 = { -4: -mp.mpf(7)/240, -3: mp.mpf(3)/10, -2: -mp.mpf(169)/120, -1: mp.mpf(61)/30,
       1: -mp.mpf(61)/30, 2: mp.mpf(169)/120, 3: -mp.mpf(3)/10, 4: mp.mpf(7)/240 }
# Richardson over several h
res={}
for hs in ['0.02','0.01','0.005','0.0025']:
    h=mp.mpf(hs)
    vals={k:perim(1+k*h,m) for k in range(-4,5)}
    res[hs]=sum(c3[k]*vals[k] for k in c3)/h**3
    print(f"12gon h={hs}: P3={mp.nstr(res[hs],35)}",flush=True)
hsl=[mp.mpf(x) for x in ['0.02','0.01','0.005','0.0025']]
M=mp.matrix([[h**0,h**2,h**4,h**6] for h in hsl]); b=mp.matrix([res[x] for x in ['0.02','0.01','0.005','0.0025']])
sol=mp.lu_solve(M,b)
print("12gon Richardson P3 =", mp.nstr(sol[0],35))
known = mp.mpf('-5.20324042951766621029582970427243461666280533')
print("12gon known      P3 =", mp.nstr(known,35))
print("diff:", mp.nstr(abs(sol[0]-known),4))
open('progress.log','a').write('PG6 s_A5: 12-gon control Richardson P3 = %s vs known -5.2032404295176662, diff %s\n'%(mp.nstr(sol[0],25),mp.nstr(abs(sol[0]-known),3)))
