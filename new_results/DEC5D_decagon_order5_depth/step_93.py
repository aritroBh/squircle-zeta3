
import pickle
import mpmath as mp
import flint
from flint import fmpz_mat

def build_basis(dps):
    mp.mp.dps = dps
    r5 = mp.sqrt(5); phi = (1+r5)/2
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
    return names, [basis[n] for n in names], r5

def lll_relation(xs, dps, scale_exp):
    N = len(xs)
    SCALE = mp.mpf(10)**scale_exp
    xi = [int(mp.nint(x*SCALE)) for x in xs]
    M = fmpz_mat(N, N+1)
    for i in range(N):
        M[i, i] = 1
        M[i, N] = xi[i]
    red = M.lll()
    best = None
    for r in range(N):
        row = [int(red[r, j]) for j in range(N+1)]
        if row[0] == 0: continue
        last = abs(row[N]); maxc = max(abs(c) for c in row[:N])
        key = (last, maxc)
        if best is None or key < best[0]:
            best = (key, row)
    return best

W5L240 = mp.mpf(pickle.load(open("w5L_250.pkl","rb"))['W5L'])
W5L420 = mp.mpf(pickle.load(open("w5L_420.pkl","rb"))['W5L'])
cands = {}
for dps, W5L, se in ((240, W5L240, 195), (420, W5L420, 370)):
    names, B, r5 = build_basis(dps)
    xs = [W5L] + B + [r5*x for x in B]
    key, row = lll_relation(xs, dps, se)
    cands[dps] = row
    s = mp.mpf(0)
    for i in range(len(xs)): s += mp.mpf(row[i])*xs[i]
    print("dps %d: maxc = %d, last coord = %d, numeric residual = %s" % (dps, key[1], key[0], mp.nstr(abs(s),3)))
print("\ncoefficient vectors identical:", cands[240] == cands[420])
if cands[240] != cands[420]:
    d = sum(1 for a,b in zip(cands[240],cands[420]) if a!=b)
    print("differing entries:", d)
