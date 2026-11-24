import mpmath as mp, pickle, sympy as sp
mp.mp.dps = 105
def eta(s): return (1-mp.power(2,1-s))*mp.zeta(s)
# I table from quadrature (fresh, extended range)
Icache = {}
def I(m, j):
    if (m,j) not in Icache:
        Icache[(m,j)] = mp.quad(lambda s: s**m*mp.log1p(mp.e**(-s))**j, [0, 0.5, 2, 8, 32, mp.inf])
    return Icache[(m,j)]
e = {}
with open('e_rm_mine.pkl','rb') as f:
    for key, val in pickle.load(f).items():
        r, m = map(int, key.split(',')); e[(r,m)] = mp.mpf(val)
def Dtilde(n, r):
    tot = mp.mpf(0)
    for j in range(1, n-2*r):
        m = n-1-j
        tot += e[(r,m)]*I(m,j)/mp.factorial(j)
    return 2*tot
print("EVEN TOPS: D~_{2n,n-1} vs 2(2n-2)! eta(2n)")
for nn in range(2, 7):
    lhs = Dtilde(2*nn, nn-1); rhs = 2*mp.factorial(2*nn-2)*eta(2*nn)
    print(f"  n={nn}: D~_{{2n,n-1}}={mp.nstr(lhs,30)}  diff={mp.nstr(abs(lhs-rhs),3)}")
print("ODD TOPS: D~_{2n+1,n-1} vs (2n-2)![z(2n+1) + 2 sum eta(2j) z(2n+1-2j)]")
for nn in range(1, 6):
    lhs = Dtilde(2*nn+1, nn-1)
    rhs = mp.factorial(2*nn-2)*(mp.zeta(2*nn+1) + 2*mp.fsum(eta(2*j)*mp.zeta(2*nn+1-2*j) for j in range(1, nn)))
    print(f"  n={nn}: D~={mp.nstr(lhs,30)}  diff={mp.nstr(abs(lhs-rhs),3)}")
print("nonzero (degree attained): all diffs tiny & values", 
      [mp.nstr(Dtilde(2*nn,nn-1),4) for nn in range(2,7)], [mp.nstr(Dtilde(2*nn+1,nn-1),4) for nn in range(1,6)])
pickle.dump(Icache, open('Icache.pkl','wb'))
