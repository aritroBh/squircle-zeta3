
# step_06.py : d7 PSLQ vs weight-7 basis; individual I_ab identification; d9; verify d5.
import mpmath as mp
mp.mp.dps = 105
import sympy as sp

Kmax = 60
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

def int_u_log(s, m, delta):
    ld = mp.log(delta); tot = mp.mpf(0)
    for i in range(m+1):
        tot += (-1)**(m-i)*mp.factorial(m)/mp.factorial(i)*ld**i/(s+1)**(m-i+1)
    return delta**(s+1)*tot

def I_ab(a, b, delta=mp.mpf('0.5')):
    near = mp.mpf(0)
    for j in range(0, b+1):
        pj = GP[j]
        for k2, coeff in enumerate(pj):
            if coeff == 0: continue
            near += mp.binomial(b,j)*coeff*int_u_log(a+k2, b-j, delta)
    far = mp.quad(lambda uu: uu**a*mp.log(1-mp.exp(-uu))**b, [delta,1,3,8,20,mp.inf])
    return near+far

need = [(3,1),(2,2),(1,3),(0,4),(5,1),(4,2),(3,3),(2,4),(1,5),(0,6),(7,1),(6,2),(5,3),(4,4),(3,5),(2,6),(1,7),(0,8),(1,2),(0,3)]
I = {}
for ab in need:
    I[ab] = I_ab(*ab)
# verify I_12 = pi^4/180 ( = 2*zeta(3,1) )
print('I_12 - pi^4/180 =', mp.nstr(I[(1,2)]-mp.pi**4/180,3))
print('I_03 - 6*zeta(2,1,1): zeta(2,1,1)=?  I_03/6 =', mp.nstr(I[(0,3)]/6,20), ' pi^4/360 =', mp.nstr(mp.pi**4/360,20), ' err', mp.nstr(abs(I[(0,3)]/6-mp.pi**4/360),3))

d5 = -(-4*I[(3,1)]+6*I[(2,2)]-4*I[(1,3)]+I[(0,4)])/2880
d7 = (-6*I[(5,1)]+15*I[(4,2)]-20*I[(3,3)]+15*I[(2,4)]-6*I[(1,5)]+I[(0,6)])/181440
d9 = -(-8*I[(7,1)]+28*I[(6,2)]-56*I[(5,3)]+70*I[(4,4)]-56*I[(3,5)]+28*I[(2,6)]-8*I[(1,7)]+I[(0,8)])/9676800
print('d5 =', mp.nstr(d5,40))
print('60*d5 +3*z5 - z2*z3 =', mp.nstr(60*d5+3*mp.zeta(5)-mp.zeta(2)*mp.zeta(3),6))
print('d7 =', mp.nstr(d7,40))
print('d9 =', mp.nstr(d9,40))

def pslq(val, basis, names, maxcoeff):
    rel = mp.pslq([val]+basis, tol=mp.mpf('1e-42'), maxcoeff=maxcoeff, maxsteps=10000)
    print(f'  PSLQ({names}, M={maxcoeff}):', rel)
    return rel

print('d7 vs {z7,z2z5,z3z4,z2^2 z3}:')
pslq(d7,[mp.zeta(7),mp.zeta(2)*mp.zeta(5),mp.zeta(3)*mp.zeta(4),mp.zeta(2)**2*mp.zeta(3)],['1','z7','z2z5','z3z4','z2^2z3'],10**8)
print('d9 vs {z9,z2z7,z3z6,z4z5,z2^2z5,z2z3z4,z3^3,z2^3z3}:')
pslq(d9,[mp.zeta(9),mp.zeta(2)*mp.zeta(7),mp.zeta(3)*mp.zeta(6),mp.zeta(4)*mp.zeta(5),mp.zeta(2)**2*mp.zeta(5),mp.zeta(2)*mp.zeta(3)*mp.zeta(4),mp.zeta(3)**3,mp.zeta(2)**3*mp.zeta(3)],['1','z9','z2z7','z3z6','z4z5','z2^2z5','z2z3z4','z3^3','z2^3z3'],10**7)

# individual I_ab identifications (weight 5 and 7 pieces)
print('individual I_ab PSLQ:')
for (a,b),basis,names in [((2,2),[mp.zeta(5),mp.zeta(2)*mp.zeta(3)],['z5','z2z3']),
                           ((1,3),[mp.zeta(5),mp.zeta(2)*mp.zeta(3)],['z5','z2z3']),
                           ((0,4),[mp.zeta(5),mp.zeta(2)*mp.zeta(3)],['z5','z2z3']),
                           ((4,2),[mp.zeta(7),mp.zeta(2)*mp.zeta(5),mp.zeta(3)*mp.zeta(4),mp.zeta(2)**2*mp.zeta(3)],['z7','z2z5','z3z4','z2^2z3']),
                           ((3,3),[mp.zeta(7),mp.zeta(2)*mp.zeta(5),mp.zeta(3)*mp.zeta(4),mp.zeta(2)**2*mp.zeta(3)],['z7','z2z5','z3z4','z2^2z3']),
                           ((2,4),[mp.zeta(7),mp.zeta(2)*mp.zeta(5),mp.zeta(3)*mp.zeta(4),mp.zeta(2)**2*mp.zeta(3)],['z7','z2z5','z3z4','z2^2z3']),
                           ((1,5),[mp.zeta(7),mp.zeta(2)*mp.zeta(5),mp.zeta(3)*mp.zeta(4),mp.zeta(2)**2*mp.zeta(3)],['z7','z2z5','z3z4','z2^2z3']),
                           ((0,6),[mp.zeta(7),mp.zeta(2)*mp.zeta(5),mp.zeta(3)*mp.zeta(4),mp.zeta(2)**2*mp.zeta(3)],['z7','z2z5','z3z4','z2^2z3'])]:
    print(f' I_{a},{b}:', mp.pslq([I[(a,b)]]+basis, tol=mp.mpf('1e-42'), maxcoeff=10**8, maxsteps=10000), names)
