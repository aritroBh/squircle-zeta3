"""Finite differences of the direct perimeter (no expansion under the integral sign)."""
import mpmath as mp, json
from jet import perimeter_direct
from polys import polys
mp.mp.dps = 70
J = json.load(open('jet_dps90_h128.json'))
K = 9; delta = mp.mpf(1)/200
for name, U in polys().items():
    ks = list(range(-K, K+1))
    vals = [perimeter_direct(U, 1+k*delta) for k in ks]
    # exact polynomial interpolation of degree 2K in d, read off Taylor coefficients
    M = mp.matrix([[ (k*delta)**j for j in range(2*K+1)] for k in ks])
    c = mp.lu_solve(M, mp.matrix(vals))
    D = [c[0], c[1], 2*c[2], 6*c[3]]
    jet = [mp.mpf(x) for x in J[name]]
    print(name, ' |FD - jet| for P,P1,P2,P3:', ['%.1e' % float(abs(D[i]-jet[i])) for i in range(4)], flush=True)
