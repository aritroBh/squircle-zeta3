
"""PG5 step_06b: gold P3(1) via roots-of-unity filter on P(p) (Taylor coefficient c3).
P(p) analytic near p=1; c3 = (1/N) sum_j zeta^{-3j} P(1+rho zeta^j)/rho^3; P3(1) = 6 c3.
"""
import mpmath as mp, pickle
mp.mp.dps = 110
m = 7
def fsum(th, p):
    return sum(mp.exp(p*mp.log(abs(mp.cos(th - i*mp.pi/m)))) for i in range(m))
def r_and_rp(th, p):
    f = fsum(th,p)
    fp = sum(-mp.sin(th - i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    r = f**(-1/p)
    rp = -(1/p)*f**(-1/p-1)*fp
    return r, rp
def perim(p):
    g = lambda th: (lambda rr: mp.sqrt(rr[0]**2+rr[1]**2))(r_and_rp(th,p))
    return 2*m*mp.quad(g, [mp.pi/2+mp.mpf('1e-9'), mp.pi/2+mp.pi/m-mp.mpf('1e-9')])
N = 16; rho = mp.mpf('0.02')
tot = mp.mpc(0)
for j in range(N):
    z = mp.exp(2*mp.pi*1j*j/N)
    Pv = perim(1+rho*z)
    tot += z**(-3)*Pv
    print(j, mp.nstr(Pv,15), flush=True)
c3 = tot/N/rho**3
P3 = 6*c3
print("c3 =", mp.nstr(c3,30))
print("P3(1) gold =", mp.nstr(P3, 40))
# repeat with rho=0.01 for error estimate
tot2 = mp.mpc(0)
for j in range(N):
    z = mp.exp(2*mp.pi*1j*j/N)
    tot2 += z**(-3)*perim(1+mp.mpf('0.01')*z)
c3b = tot2/N/mp.mpf('0.01')**3
print("P3 (rho=0.01) =", mp.nstr(6*c3b,40))
print("agreement:", mp.nstr(abs(P3-6*c3b),5))
pickle.dump({'P3_gold':str(P3),'P3_rho2':str(6*c3b)}, open('gold14.pkl','wb'))
open('progress.log','a').write('PG5 step_06b: gold P3(1) 14-gon = %s vs rho2 %s\n'%(mp.nstr(P3,30),mp.nstr(abs(P3-6*c3b),3)))
