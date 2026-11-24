
# step_13.py : granular central-binomial building-block sums for the triangle end; PSLQ each.
import mpmath as mp
mp.mp.dps = 100

def f(m): return mp.gamma(m+1)**2/mp.gamma(2*m+1)
N = 800
def H(m,j): return mp.nsum(lambda k: k**(-j), [1,m]) if m>0 else mp.mpf(0)
# faster: cumulative
maxm = N
H1=[mp.mpf(0)]; H2=[mp.mpf(0)]; H3=[mp.mpf(0)]
for m in range(1,maxm+1):
    H1.append(H1[-1]+mp.mpf(1)/m); H2.append(H2[-1]+mp.mpf(1)/m**2); H3.append(H3[-1]+mp.mpf(1)/m**3)
sums = {}
def acc(name, term):
    sums.setdefault(name, mp.mpf(0))
    sums[name] += term
for m in range(1, N+1):
    fm = f(m)
    if fm < mp.mpf('1e-115') and m > 60: break
    acc('S0', fm)                       # sum f(m)
    acc('S1', fm/m)                     # = zeta2/3
    acc('Sm', m*fm)                     # sum m f(m)
    acc('T1', fm*(H1[2*m]-H1[m])/m)     # c1 = 2 T1
    acc('P2', fm*H2[m])
    acc('P3', fm*H2[2*m])
    acc('P4', fm*(H1[2*m]-H1[m])**2)
    acc('Q2', m*fm*H3[m])
    acc('Q3', m*fm*H3[2*m])
    acc('Q4', m*fm*H1[m])               # weight-4? m*fm*H_m
rt3 = mp.sqrt(3); Cl = mp.clsin(2, mp.pi/3)
for nm in sums:
    print(nm, '=', mp.nstr(sums[nm], 30))
basis = [mp.mpf(1), mp.pi/rt3, mp.log(2), mp.log(3), mp.zeta(2), mp.zeta(3), Cl, Cl*mp.pi]
bn = ['1','pi/rt3','log2','log3','z2','z3','Cl','pi*Cl']
for nm in sums:
    rel = mp.pslq([sums[nm]]+basis, tol=mp.mpf('1e-38'), maxcoeff=10**6, maxsteps=10000)
    print(f'{nm} PSLQ:', rel)
