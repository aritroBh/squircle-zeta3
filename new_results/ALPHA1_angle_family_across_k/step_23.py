"""step_23.py -- weight-three slope value: slope(alpha) = s(alpha) J(r), identify J / slope in c."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 80

def Jval(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def lmt(w): return mp.log((1-w)*(w+r)/(2*(1-r)*w))
    def f(w):
        return (1/(w+1) - 1/(w+r))*lt(w)**2 + (1/w - 2/(w+1))*lt(w)*lmt(w)
    return mp.quad(f, [r, mp.sqrt(r), 1])

def slope(c):
    return (1-c**2)/(2*c**2) * Jval(c)

# squircle control
lam = mp.log(2)
sq = -(mp.zeta(3)/8 + mp.pi**2*lam/24 + lam**3/12)
print("squircle slope =", mp.nstr(slope(1/mp.sqrt(2)), 30), " target:", mp.nstr(sq, 30),
      " diff:", mp.nstr(abs(slope(1/mp.sqrt(2))-sq), 3))

angles = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7, mp.mpf(5)/6, mp.mpf(1)/5, mp.mpf(3)/7]
SL = [slope(c) for c in angles]
for c, v in zip(angles, SL):
    print(f"c={mp.nstr(c,8)}: slope={mp.nstr(v,35)}")

# identify: basis of weight-3 constants in c: zeta3, pi^2 logs, log^3, Li3 at c, -c, (1-c)/2, ...
def basis(c):
    L = mp.log; A = lambda x: mp.polylog(3, x)
    return [('zeta3', mp.zeta(3)),
            ('pi^2 log2', mp.pi**2*L(2)), ('pi^2 logc', mp.pi**2*L(c)), ('pi^2 log(1+c)', mp.pi**2*L(1+c)), ('pi^2 log(1-c)', mp.pi**2*L(1-c)),
            ('log2^3', L(2)**3), ('logc^3', L(c)**3), ('log(1+c)^3', L(1+c)**3), ('log(1-c)^3', L(1-c)**3),
            ('log2^2 logc', L(2)**2*L(c)), ('log2 logc^2', L(2)*L(c)**2), ('logc^2 log(1+c)', L(c)**2*L(1+c)),
            ('logc log(1+c)^2', L(c)*L(1+c)**2), ('log(1+c)^2 log(1-c)', L(1+c)**2*L(1-c)),
            ('Li3(c)', A(c)), ('Li3(-c)', A(-c)), ('Li3((1-c)/2)', A((1-c)/2)), ('Li3((1+c)/2)', A((1+c)/2)),
            ('Li3(1-c)', A(1-c)), ('Li3(2c/(1+c))', A(2*c/(1+c))), ('Li3(c/(1+c))', A(c/(1+c)))]

C = [[x for _, x in basis(c)] for c in angles]
names = [n for n, _ in basis(angles[0])]
M, N = len(angles), len(names)
def n2(v): return mp.sqrt(sum(x*x for x in v))
A = [mp.matrix([C[i][j] for i in range(M)]) for j in range(N)]
for j in range(N): A[j] = A[j]/n2(A[j])
sel = []
for j in range(N):
    v = A[j].copy()
    for s in sel: v = v - (v.T*A[s])[0]*A[s]
    if n2(v) > mp.mpf('1e-50'):
        for s in sel: v = v - (v.T*A[s])[0]*A[s]
        v = v/n2(v); A[j] = v; sel.append(j)
Asel = mp.matrix(M, len(sel))
for k, j in enumerate(sel):
    for i in range(M): Asel[i, k] = C[i][j]
Tv = mp.matrix(SL)
q = mp.lu_solve(Asel.T*Asel, Asel.T*Tv)
resid = Asel*q - Tv
print("LS residual:", mp.nstr(max(abs(x) for x in resid), 3))
coef = {}
for k, j in enumerate(sel):
    if abs(q[k]) > mp.mpf('1e-50'):
        coef[names[j]] = q[k]
        print(f"   {names[j]}: {mp.nstr(q[k], 35)}")
for c in (mp.mpf(1)/7, mp.mpf(4)/9, 1/mp.sqrt(2)):
    bmap = {n: x for n, x in basis(c)}
    pred = sum(coef[n]*bmap[n] for n in coef)
    print(f"   holdout c={mp.nstr(c,9)}: diff={mp.nstr(abs(slope(c)-pred),3)}")
print("STEP23 DONE")
