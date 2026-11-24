
# step_07.py : corrected PSLQ bases; exact d7, d9; individual I_ab table; MZV sanity checks.
import mpmath as mp
mp.mp.dps = 105
import sympy as sp

Kmax = 70
u = sp.symbols('u')
ex = sp.log(sum((-u)**j/sp.factorial(j+1) for j in range(Kmax+3)))
cs = [sp.expand(ex).series(u,0,Kmax+1).removeO().coeff(u,k) for k in range(Kmax+1)]
c = [mp.mpf(x) for x in cs]
def gpow(j):
    poly = [mp.mpf(1)]
    for _ in range(j):
        new = [mp.mpf(0)]*(len(poly)+Kmax)
        for i,a in enumerate(poly):
            for k2 in range(1, Kmax+1): new[i+k2] += a*c[k2]
        poly = new
    return poly
GP = {j: gpow(j) for j in range(0,11)}
def int_u_log(s,m,delta):
    ld = mp.log(delta); tot = mp.mpf(0)
    for i in range(m+1):
        tot += (-1)**(m-i)*mp.factorial(m)/mp.factorial(i)*ld**i/(s+1)**(m-i+1)
    return delta**(s+1)*tot
def I_ab(a,b,delta=mp.mpf('0.5')):
    near = mp.mpf(0)
    for j in range(0,b+1):
        for k2,coeff in enumerate(GP[j]):
            if coeff==0: continue
            near += mp.binomial(b,j)*coeff*int_u_log(a+k2,b-j,delta)
    far = mp.quad(lambda uu: uu**a*mp.log(1-mp.exp(-uu))**b,[delta,1,3,8,20,mp.inf])
    return near+far
need7 = [(5,1),(4,2),(3,3),(2,4),(1,5),(0,6)]
need9 = [(7,1),(6,2),(5,3),(4,4),(3,5),(2,6),(1,7),(0,8)]
I = {ab: I_ab(*ab) for ab in need7+need9}

# sanity: direct MZV values
print('I_03/(-6) =', mp.nstr(-I_ab(0,3)/6, 15), ' vs pi^4/90 =', mp.nstr(mp.pi**4/90,15))

W7 = [mp.zeta(7), mp.zeta(2)*mp.zeta(5), mp.zeta(2)**2*mp.zeta(3)]
names7 = ['z7','z2z5','z2^2z3']
print('--- weight-7 individual I_ab:')
rel7 = {}
for ab in need7:
    r = mp.pslq([I[ab]]+W7, tol=mp.mpf('1e-42'), maxcoeff=10**8, maxsteps=10000)
    rel7[ab] = r
    print(f' I_{ab}:', r)
d7 = (-6*I[(5,1)]+15*I[(4,2)]-20*I[(3,3)]+15*I[(2,4)]-6*I[(1,5)]+I[(0,6)])/181440
print('d7 =', mp.nstr(d7,40))
print('d7 PSLQ:', mp.pslq([d7]+W7, tol=mp.mpf('1e-42'), maxcoeff=10**8, maxsteps=10000), names7)

# weight 9: products and their relation
z9, z2z7, p6z3, p4z5, p4z2z3, z33 = mp.zeta(9), mp.zeta(2)*mp.zeta(7), mp.pi**6*mp.zeta(3), mp.pi**4*mp.zeta(5), mp.pi**4*mp.zeta(2)*mp.zeta(3), mp.zeta(3)**3
prods = [z9, z2z7, p6z3, p4z5, p4z2z3, z33]
pnames = ['z9','z2z7','pi^6z3','pi^4z5','pi^4z2z3','z3^3']
print('relation among weight-9 products:', mp.pslq(prods, tol=mp.mpf('1e-42'), maxcoeff=10**8, maxsteps=10000), pnames)
d9 = -(-8*I[(7,1)]+28*I[(6,2)]-56*I[(5,3)]+70*I[(4,4)]-56*I[(3,5)]+28*I[(2,6)]-8*I[(1,7)]+I[(0,8)])/9676800
print('d9 =', mp.nstr(d9,40))
print('d9 PSLQ vs products:', mp.pslq([d9]+prods, tol=mp.mpf('1e-42'), maxcoeff=10**8, maxsteps=10000), ['1']+pnames)
