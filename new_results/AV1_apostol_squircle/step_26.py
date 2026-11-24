code = r'''
# step_22.py : d11 and weight-11 individuals vs the 7 independent products.
import mpmath as mp
mp.mp.dps = 130
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
need=[(10-j,j) for j in range(1,11)]
I={ab:I_ab(*ab) for ab in need}
B10=mp.mpf(5)/66
d11 = B10/(10*mp.factorial(10))*sum(mp.binomial(10,j)*(-1)**j*I[(10-j,j)] for j in range(1,11))
z=mp.zeta
W11=[z(11),z(2)*z(9),mp.pi**4*z(7),mp.pi**6*z(5),mp.pi**8*z(3),z(3)**2*z(5),z(2)*z(3)**3]
n11=['z11','z2z9','pi4z7','pi6z5','pi8z3','z3^2z5','z2z3^3']
print('individuals:')
for ab in need:
    r = mp.pslq([I[ab]]+W11,tol=mp.mpf('1e-60'),maxcoeff=10**10,maxsteps=10000)
    print(f' I_{ab}:', r)
rel = mp.pslq([d11]+W11,tol=mp.mpf('1e-60'),maxcoeff=10**10,maxsteps=10000)
print('d11 PSLQ:', rel)
if rel and rel[0]!=0:
    print(' => d11 = -(', '+'.join(f'{rel[i+1]}*{n11[i]}' for i in range(7)), f')/{rel[0]}')
    val=-sum(mp.mpf(rel[i+1])/rel[0]*b for i,b in enumerate(W11))
    print(' residual:', mp.nstr(abs(d11-val),4))
'''
open('steps/step_22.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_22.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-800:] if r.returncode else '')
open('progress.log','a').write('step_22: d11 identified vs 7 independent weight-11 products\n')
