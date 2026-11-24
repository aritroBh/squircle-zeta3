
import pickle
import mpmath as mp
mp.mp.dps = 240
import flint
from flint import fmpz, fmpz_mat
r5 = mp.sqrt(5); phi = (1+r5)/2
W5L = mp.mpf(pickle.load(open("w5L_250.pkl","rb"))['W5L'])
l2 = mp.log(2); lf = mp.log(phi)
z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
L3p = mp.polylog(3, 1/phi); L3m = mp.polylog(3, -1/phi)
L4p = mp.polylog(4, 1/phi); L4m = mp.polylog(4, -1/phi)
L5p = mp.polylog(5, 1/phi); L5m = mp.polylog(5, -1/phi)
basis = {'z5': z5, 'z2z3': z2*z3, 'L5p': L5p, 'L5m': L5m,
         'L4p l2': L4p*l2, 'L4p lf': L4p*lf, 'L4m l2': L4m*l2, 'L4m lf': L4m*lf,
         'L3p l2^2': L3p*l2**2, 'L3p l2lf': L3p*l2*lf, 'L3p lf^2': L3p*lf**2,
         'L3m l2^2': L3m*l2**2, 'L3m l2lf': L3m*l2*lf, 'L3m lf^2': L3m*lf**2,
         'z3 l2^2': z3*l2**2, 'z3 l2lf': z3*l2*lf, 'z3 lf^2': z3*lf**2,
         'z2^2 l2': z2**2*l2, 'z2^2 lf': z2**2*lf,
         'z2 l2^3': z2*l2**3, 'z2 l2^2lf': z2*l2**2*lf, 'z2 l2lf^2': z2*l2*lf**2, 'z2 lf^3': z2*lf**3,
         'l2^5': l2**5, 'l2^4lf': l2**4*lf, 'l2^3lf^2': l2**3*lf**2,
         'l2^2lf^3': l2**2*lf**3, 'l2lf^4': l2*lf**4, 'lf^5': lf**5, '1': mp.mpf(1)}
names = list(basis)
B = [basis[n] for n in names]
xs = [W5L] + B + [r5*x for x in B]   # find relation with W5L = sum (p_i + q_i sqrt5) B_i
# but we want p,q separated: unknowns c_i (for B_i) and d_i (for sqrt5 B_i); target row 0 = W5L.
# Standard integer relation among all xs: find m: sum m_i x_i = 0 with m_0 != 0.
N = len(xs)
SCALE = 10**195
xi = [int(mp.nint(x*SCALE)) for x in xs]
# rows: e_i + xi[i] e_{N+1}
M = fmpz_mat(N, N+1)
for i in range(N):
    M[i, i] = 1
    M[i, N] = xi[i]
import time
t0 = time.time()
red = M.lll()
print("LLL in %.1fs" % (time.time()-t0))
row = [int(red[0, j]) for j in range(N)]
last = int(red[0, N])
print("first reduced row head:", row[:6], "... max|c| =", max(abs(c) for c in row))
print("last coord:", last)
rel = row
if rel[0] != 0:
    s = sum(mp.mpf(rel[i])*xs[i] for i in range(N))
    print("numeric residual of relation:", mp.nstr(abs(s), 3))
    # express W5L
    p = [-mp.mpf(rel[i])/mp.mpf(rel[0]) for i in range(1, N)]
    s2 = sum(p[i-1]*xs[i] for i in range(1, N))
    print("W5L vs combination:", mp.nstr(abs(s2 - W5L), 3))
    print("p-part:")
    for i, n in enumerate(names):
        if rel[i+1]: print("   %-12s %s/(%s)" % (n, -rel[i+1], rel[0]))
    print("q-part:")
    for i, n in enumerate(names):
        j = rel[len(names)+1+i]
        if j: print("   %-12s %s/(%s)" % (n, -j, rel[0]))
