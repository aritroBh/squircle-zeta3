code = r'''
# step_09.py : resolve d9 — exact rational assembly from individual I_ab relations vs direct PSLQ.
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = 105
import sympy as sp
Kmax = 70
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
GP={j:gpow(j) for j in range(0,11)}
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
    return near+mp.quad(lambda uu: uu**a*mp.log(1-mp.exp(-uu))**b,[delta,1,3,8,20,mp.inf])
need9=[(7,1),(6,2),(5,3),(4,4),(3,5),(2,6),(1,7),(0,8)]
I={ab:I_ab(*ab) for ab in need9}
# individual relations as exact rationals: I_ab = sum r_j * basis_j, basis = [z9,z2z7,pi^6z3,pi^4z5,z3^3]
rels = {(7,1):[-5040,0,0,0,0],
        (6,2):[Fr(120960,21),Fr(-30240,21),Fr(-32,21),Fr(-336,21),0],
        (5,3):[Fr(-20160,3),Fr(6480,3),Fr(4,3),Fr(54,3),Fr(-360,3)],
        (4,4):[Fr(-120960,15),Fr(43200,15),Fr(20,15),Fr(336,15),Fr(-4320,15)],
        (3,5):[Fr(-30240,3),Fr(10800,3),Fr(5,3),Fr(84,3),Fr(-1080,3)],
        (2,6):[Fr(-40320,3),Fr(12960,3),Fr(8,3),Fr(108,3),Fr(-720,3)],
        (1,7):[Fr(-60480,3),Fr(15120,3),Fr(16,3),Fr(168,3),0],
        (0,8):[40320,0,0,0,0]}
rels = {k:[Fr(x) for x in v] for k,v in rels.items()}
# verify each relation numerically
zb = [mp.zeta(9),mp.zeta(2)*mp.zeta(7),mp.pi**6*mp.zeta(3),mp.pi**4*mp.zeta(5),mp.zeta(3)**3]
for ab,r in rels.items():
    val = sum(mp.mpf(x)*b for x,b in zip(r,zb))
    print(f'I_{ab} relation residual: {mp.nstr(abs(I[ab]-val),3)}')
# assemble d9 = -( -8I71+28I62-56I53+70I44-56I35+28I26-8I17+I08 )/9676800
comb = [( (7,1),-8),((6,2),28),((5,3),-56),((4,4),70),((3,5),-56),((2,6),28),((1,7),-8),((0,8),1)]
d9r = [Fr(0)]*5
for ab,w in comb:
    for j in range(5): d9r[j] -= w*rels[ab][j]
d9r = [x/Fr(9676800) for x in d9r]
print('d9 assembled coeffs [z9,z2z7,pi^6z3,pi^4z5,z3^3]:', [str(x) for x in d9r])
d9asm = sum(mp.mpf(x)*b for x,b in zip(d9r,zb))
d9num = -(-8*I[(7,1)]+28*I[(6,2)]-56*I[(5,3)]+70*I[(4,4)]-56*I[(3,5)]+28*I[(2,6)]-8*I[(1,7)]+I[(0,8)])/9676800
print('d9 numeric     =', mp.nstr(d9num,40))
print('d9 assembled   =', mp.nstr(d9asm,40))
print('d9 direct PSLQ =', mp.nstr((-428400*zb[0]+136080*zb[1]+79*zb[2]+1134*zb[3]-10080*zb[4])/1814400,40))
print('|num - assembled| =', mp.nstr(abs(d9num-d9asm),3))
'''
open('steps/step_09.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_09.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')
open('progress.log','a').write('step_09: d9 assembly cross-check\n')
