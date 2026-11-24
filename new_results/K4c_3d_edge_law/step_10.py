"""step_10.py -- PSLQ implementation; analyze q = S1/S0 for models A and B."""
import os, json
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf

def pslq(x, tol=None, maxiter=10000):
    """classic PSLQ. x: list of mpf. Returns integer relation vector a with sum a_i x_i ~ 0."""
    n = len(x)
    dps = mp.mp.dps
    if tol is None:
        tol = mpf(10) ** (-(dps * 3 // 4))
    # scale
    A = mp.matrix(n)
    for i in range(n):
        A[i,i] = 1
    # H matrix (n x n-1)
    H = mp.matrix(n, n-1)
    for j in range(n-1):
        s = mp.sqrt(sum(x[i]**2 for i in range(j, n)))
        if x[j] < 0: s = -s
        H[j,j] = s
        for i in range(j+1, n):
            H[i,j] = x[i] * x[j] / (s * (s + abs(x[j]))) if False else 0
        # standard: use the standard PSLQ initialization
        for i in range(j+1, n):
            H[i,j] = x[i] * (x[j] / s) / (s + abs(x[j]))
        H[j,j] = - (x[j]/s) * s / (s + abs(x[j])) if False else H[j,j]
    # Actually implement standard: re-derive below
    return None

# Use the well-tested formulation instead:
def pslq2(x, tol=None, maxcoeff=None, maxsteps=100000):
    """PSLQ following Bailey-Broadhurst. Returns relation list or None."""
    n = len(x)
    gam = mp.sqrt(mpf(4)/3)
    # y = x / norm(x)
    nx = mp.sqrt(sum(xi**2 for xi in x))
    y = [xi/nx for xi in x]
    # A = I, B = I
    A = mp.eye(n)
    B = mp.eye(n)
    # H: lower trapezoidal (n x n-1) via partial LQ decomposition
    H = mp.matrix(n, n-1)
    for i in range(n):
        for j in range(min(i+1, n-1)):
            H[i,j] = 0
    for j in range(n-1):
        # compute s_j = sqrt(sum_{i>=j} y_i^2) using current y (after previous steps)
        pass
    # build H directly from y:
    yv = y[:]
    for j in range(n-1):
        s = mp.sqrt(sum(yv[i]**2 for i in range(j, n)))
        if yv[j] < 0: s = -s
        # v vector
        v = [mpf(0)]*n
        v[j] = yv[j] + s
        for i in range(j+1, n): v[i] = yv[i]
        # apply Householder to yv
        vv = sum(t**2 for t in v)
        if vv == 0: continue
        fac = 2*sum(v[i]*yv[i] for i in range(n))/vv
        for i in range(n): yv[i] -= fac*v[i]
        # set H column j: H[i,j] = (householder applied)
        for i in range(n):
            pass
    # simpler: compute H by QR of column vector y (n x 1): H = R of [y]? Standard PSLQ:
    # Instead implement the H update directly:
    H = mp.matrix(n, n-1)
    s = [mpf(0)]*n
    for k in range(n-1, -1, -1):
        s[k] = mp.sqrt(sum(y[i]**2 for i in range(k, n)))
    for j in range(n-1):
        for i in range(n):
            if i < j: H[i,j] = 0
            elif i == j:
                H[i,j] = s[j+1]/s[j]
            else:
                H[i,j] = -y[i]*y[j]/(s[j]*s[j+1])
    # main iteration
    if tol is None:
        tol = mpf(10)**(-mp.mp.dps*2/3)
    for step in range(maxsteps):
        # find r maximizing gam^r |H[r,r]|
        best = -1; r = 0
        for i in range(n-1):
            val = gam**i * abs(H[i,i])
            if val > best: best = val; r = i
        # swap rows r, r+1 of H, A, B... and y
        H.swap_row(r, r+1)
        # corner restore
        if r < n-2:
            t0 = mp.sqrt(H[r,r]**2 + H[r,r+1]**2)
            t1 = H[r,r]/t0; t2 = H[r,r+1]/t0
            for i in range(r, n):
                t3 = H[i,r]; t4 = H[i,r+1]
                H[i,r] = t1*t3 + t2*t4
                H[i,r+1] = -t2*t3 + t1*t4
        # update y, A, B
        y[r], y[r+1] = y[r+1], y[r]
        A.swap_row(r, r+1)
        B.swap_col(r, r+1)
        # reduce: for i > r+1
        for i in range(r+1, n):
            for j in range(min(i, n-2), -1, -1):
                t = mp.nint(H[i,j]/H[j,j]) if H[j,j] != 0 else 0
                if t != 0:
                    y[j] += t*y[i]
                    for k in range(j+1):
                        H[i,k] -= t*H[j,k]
                    for k in range(n):
                        A[i,k] -= t*A[j,k]
                        B[k,j] += t*B[k,i]
        # check termination: norm of H column? Standard: relation found if some y_i ~ 0
        for i in range(n):
            if abs(y[i]) < tol:
                # relation is column i of B
                rel = [B[j,i] for j in range(n)]
                return [int(mp.re(t)) for t in rel], abs(y[i])
    return None, None

if __name__ == "__main__":
    mp.mp.dps = 60
    l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
    data = {}
    for name, fn in [("A","modelA_h96.json"), ("B","modelB_h48.json")]:
        J = json.load(open(fn)); S = [mpf(s) for s in J["Sdk"]]
        data[name] = S
        q = S[1]/S[0]
        print(f"model {name}: q = {mp.nstr(q,45)}")
        for basis_name, basis in [("1,l2,l3", [mpf(1), l2, l3]),
                                  ("1,l2,l3,l5", [mpf(1), l2, l3, l5]),
                                  ("l2,l3", [l2, l3]),
                                  ("1,l2,l3,l5,l7", [mpf(1), l2, l3, l5, mp.log(7)])]:
            rel, res = pslq2([q] + basis, tol=mpf('1e-40'))
            if rel:
                # q = -(sum rel[1:] . basis)/rel[0]
                from fractions import Fraction
                terms = [Fraction(rel[i+1]).limit_denominator(10**9) for i in range(len(basis))]
                print(f"   basis[{basis_name}]: rel={rel} resid={mp.nstr(res,3)}")
            else:
                print(f"   basis[{basis_name}]: none found")
