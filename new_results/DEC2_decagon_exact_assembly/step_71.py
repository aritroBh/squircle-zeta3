import mpmath as mp
mp.mp.dps = 80
z2 = mp.pi**2/6; z3 = mp.zeta(3)
Li3 = lambda x: mp.polylog(3,x)
def basis(x):
    lx = mp.log(x); l1x = mp.log(1-x)
    return [l1x**3, l1x**2*lx, l1x*lx**2, lx**3, l1x**2, l1x*lx, lx**2]
xs = [mp.mpf(v)/100 for v in [17, 23, 31, 41, 53, 61, 71]]
A = []
b = []
for x in xs:
    LHS = Li3(x)+Li3(1-x)+Li3(x/(x-1))
    base = LHS - (z3 + z2*mp.log(1-x))
    A.append(basis(x))
    b.append(base)
# solve A c = b
M = mp.matrix(A); bb = mp.matrix(b)
c = mp.lu_solve(M, bb)
for row in c:
    print(mp.nstr(row, 30), " ~ ", mp.nstr(mp.pslq([row, mp.mpf(1), mp.mpf(2), mp.mpf(3), mp.mpf(6)], maxcoeff=1000000) if False else row, 12))
