import mpmath as mp
mp.mp.dps=80
z3=mp.zeta(3);pi=mp.pi
I=mp.quad(lambda th: mp.log(mp.sin(th))**3,[0,mp.pi/2],method='tanh-sinh')
# exact known: -(pi/2)log^3 2 - (pi^3/8)log2 - (3pi/8)?? compute residual
for basis,bn in [([z3,pi**3,pi**2*mp.log(2),pi*mp.log(2)**2,mp.log(2)**3],['z3','pi^3','pi^2log2','pi*log2^2','log2^3'])]:
    rel=mp.pslq(mp.matrix([I]+basis),tol=mp.mpf('1e-75'),maxcoeff=1000000)
    print('plain pslq:',rel)
    if rel is not None and rel[0]!=0:
        print('  I =', mp.nstr(-sum(mp.mpf(rel[k+1])*basis[k] for k in range(len(basis)))/rel[0],30))
        print('  z3 coeff=',mp.nstr(-mp.mpf(rel[1])/rel[0],8))
# residual after removing -(pi/2)log^3 2 - (pi^3/8)log2:
r1 = I + (pi/2)*mp.log(2)**3 + (pi**3/8)*mp.log(2)
print('residual r1=',mp.nstr(r1,25),' r1/z3=',mp.nstr(r1/z3,12),' r1/(pi*z3)=',mp.nstr(r1/(pi*z3),12))
