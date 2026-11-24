code = r'''
# step_20.py : (i) d1 PSLQ with correct 85 digits; (ii) d11 via I_ab (weight 11) + PSLQ.
import mpmath as mp
mp.mp.dps = 130
d1 = mp.mpf('1.73300259020287278417012270424573195250974925287880319940547987633644278103184087570877319')
import sympy as sp
Kmax = 80
u = sp.symbols('u')
ex = sp.log(sum((-u)**j/sp.factorial(j+1) for j in range(Kmax+3)))
cs=[sp.expand(ex).series(u,0,Kmax+1).removeO().coeff(u,k) for k in range(Kmax+1)]
c=[mp.mpf(x) for x in cs]
def gpow(j):
    poly=[mp.mpf(1)]
    for _ in range(j):
        new=[mp.mpf(0)]*(len(poly)+Kmax)
        for i,a in enumerate(poly):
            for k2 in range(1,Kmax+1): new[i+k2]+=a*c[k2]
        poly=new
    return poly
GP={j:gpow(j) for j in range(0,12)}
def int_u_log(s,m,delta):
    ld=mp.log(delta); tot=mp.mpf(0)
    for i in range(m+1): tot+=(-1)**(m-i)*mp.factorial(m)/mp.factorial(i)*ld**i/(s+1)**(m-i+1)
    return delta**(s+1)*tot
def I_ab(a,b,delta=mp.mpf('0.5')):
    near=mp.mpf(0)
    for j in range(0,b+1):
        for k2,coeff in enumerate(GP[j]):
            if coeff==0: continue
            near+=mp.binomial(b,j)*coeff*int_u_log(a+k2,b-j,delta)
    return near+mp.quad(lambda uu: uu**a*mp.log(1-mp.exp(-uu))**b,[delta,1,3,8,20,60,200])

print('--- d1 PSLQ pools (correct digits)')
gl = mp.mpf('1.2824271291006226368753425688697917277676889273250')  # Glaisher
for pool,names in [([mp.mpf(1),mp.euler,mp.log(2*mp.pi),mp.log(2),mp.zeta(2),mp.zeta(3)],'1,gam,log2pi,log2,z2,z3'),
                   ([mp.log(2),mp.log(2*mp.pi),gl,mp.euler*mp.log(2)],'log2,log2pi,Glaisher,gam*log2'),
                   ([mp.pi,mp.euler,mp.log(2),mp.zeta(2)],'pi,gam,log2,z2'),
                   ([mp.zeta(2),mp.zeta(3),mp.log(2)**2,mp.log(2)**3],'z2,z3,log2^2,log2^3')]:
    rel = mp.pslq([d1]+pool, tol=mp.mpf('1e-60'), maxcoeff=10**10, maxsteps=10000)
    print(f' d1 vs {names}:', rel)

print('--- d11')
need=[(10-j,j) for j in range(1,11)]
I={ab:I_ab(*ab) for ab in need}
import math
B10 = mp.mpf(5)/66
d11 = B10/(10*mp.factorial(10))*sum(mp.binomial(10,j)*(-1)**j*I[(10-j,j)] for j in range(1,11))
print('d11 =', mp.nstr(d11,50))
z=mp.zeta
prods=[z(11),z(2)*z(9),mp.pi**4*z(7),mp.pi**6*z(5),mp.pi**8*z(3),z(2)**2*z(7),mp.pi**4*z(2)*z(5),z(3)**2*z(5),z(2)**3*z(5),mp.pi**4*z(2)**2*z(3),z(2)*z(3)**3]
pn=['z11','z2z9','pi4z7','pi6z5','pi8z3','z2^2z7','pi4z2z5','z3^2z5','z2^3z5','pi4z2^2z3','z2z3^3']
print('relations among weight-11 products (first two found):')
print(mp.pslq(prods,tol=mp.mpf('1e-60'),maxcoeff=10**10,maxsteps=10000))
rel=mp.pslq([d11]+prods,tol=mp.mpf('1e-60'),maxcoeff=10**10,maxsteps=10000)
print('d11 PSLQ:', rel)
if rel and rel[0]!=0:
    print(' => d11 = -(', '+'.join(f'{rel[i+1]}*{pn[i]}' for i in range(11)), f')/{rel[0]}')
    q=rel[0]
    val=-sum(mp.mpf(rel[i+1])/q*b for i,b in enumerate(prods))
    print(' residual:', mp.nstr(abs(d11-val),4))
'''
open('steps/step_20.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_20.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-800:] if r.returncode else '')
open('progress.log','a').write('step_20: d1 PSLQ (no relation), d11 weight-11 PSLQ\n')
