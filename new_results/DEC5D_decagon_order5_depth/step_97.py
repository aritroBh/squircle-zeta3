
import pickle, time
import mpmath as mp
import flint
from flint import fmpz_mat
mp.mp.dps = 420
chi = {1:1, 2:-1, 3:-1, 4:1}
r5 = mp.sqrt(5); phi = (1+r5)/2
l2 = mp.log(2); lf = mp.log(phi)
z2 = mp.pi**2/6; z3 = mp.zeta(3); z5 = mp.zeta(5)
L3p = mp.polylog(3, 1/phi); L3m = mp.polylog(3, -1/phi)
L4p = mp.polylog(4, 1/phi); L4m = mp.polylog(4, -1/phi)
L5p = mp.polylog(5, 1/phi); L5m = mp.polylog(5, -1/phi)
print("atoms ready", flush=True)

def lll_best(xs, scale_exp):
    N = len(xs)
    SCALE = mp.mpf(10)**scale_exp
    xi = [int(mp.nint(x*SCALE)) for x in xs]
    M = fmpz_mat(N, N+1)
    for i in range(N):
        M[i, i] = 1
        M[i, N] = xi[i]
    red = M.lll(delta=mp.mpf('0.99'))
    best = None
    for r in range(N):
        row = [int(red[r, j]) for j in range(N+1)]
        if row[0] == 0: continue
        last = abs(row[N]); maxc = max(abs(c) for c in row[:N])
        key = (last, maxc)
        if best is None or key < best[0]: best = (key, row)
    return best

def logs_basis(L3, L4, L5, tag):
    """basis for a block with single Li atoms L3,L4,L5 (mpf) or None."""
    b = {}
    b['z5'] = z5; b['z2z3'] = z2*z3
    if L5 is not None: b['L5'+tag] = L5
    if L4 is not None:
        b['L4'+tag+' l2'] = L4*l2; b['L4'+tag+' lf'] = L4*lf
    if L3 is not None:
        b['z2 L3'+tag] = z2*L3
        b['L3'+tag+' l2^2'] = L3*l2**2; b['L3'+tag+' l2lf'] = L3*l2*lf; b['L3'+tag+' lf^2'] = L3*lf**2
    b['z3 l2^2'] = z3*l2**2; b['z3 l2lf'] = z3*l2*lf; b['z3 lf^2'] = z3*lf**2
    b['z2^2 l2'] = z2**2*l2; b['z2^2 lf'] = z2**2*lf
    b['z2 l2^3'] = z2*l2**3; b['z2 l2^2lf'] = z2*l2**2*lf; b['z2 l2lf^2'] = z2*l2*lf**2; b['z2 lf^3'] = z2*lf**3
    b['l2^5'] = l2**5; b['l2^4lf'] = l2**4*lf; b['l2^3lf^2'] = l2**3*lf**2
    b['l2^2lf^3'] = l2**2*lf**3; b['l2lf^4'] = l2*lf**4; b['lf^5'] = lf**5
    b['1'] = mp.mpf(1)
    return b

W = pickle.load(open("w5L_420.pkl","rb"))
for blk, L3, L4, L5, tag in (('single-1', L3m, L4m, L5m, 'm'),
                             ('single-4', L3p+L3m, L4p+L4m, L5p+L5m, 's')):
    target = mp.mpf(W['blocks'][blk])
    bb = logs_basis(L3, L4, L5, tag)
    names = list(bb); B = [bb[n] for n in names]
    xs = [target] + B + [r5*x for x in B]
    t0 = time.time()
    key, row = lll_best(xs, 360)
    s = sum(mp.mpf(row[i])*xs[i] for i in range(len(xs)))
    print("%s (%d elts): LLL %.1fs, maxc=%d, last=%d, resid=%s" %
          (blk, len(names), time.time()-t0, key[1], key[0], mp.nstr(abs(s),3)), flush=True)
    if key[0] < 100 and abs(s) < mp.mpf('1e-300'):
        print("  CANDIDATE FOUND")
        for i, n in enumerate(names):
            if row[i+1]: print("   %-12s (%s)/(%s)" % (n, -row[i+1], row[0]))
        for i, n in enumerate(names):
            j = row[len(names)+1+i]
            if j: print("   %-12s (%s)/(%s) *sqrt5" % (n, -j, row[0]))
