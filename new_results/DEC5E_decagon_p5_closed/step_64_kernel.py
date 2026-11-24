"""step_64_kernel.py - exact relations among the 37 invisible atoms (LLL kernel), RREF over Q.

Lattice: rows [e_i | C*v_i, C*r5*v_i] for each atom value v_i (sqrt5-paired columns),
so an integer row (m_i) with small norm gives sum m_i v_i ~ 0 AND sum m_i r5 v_i ~ 0 simultaneously,
i.e. relations with integer coefficients valid over Q(sqrt5) automatically.
Usage: python step_64_kernel.py D TAB
"""
import sys, os, pickle, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from flint import fmpz_mat
HERE = os.path.dirname(os.path.abspath(__file__))

def setup(D):
    mp.mp.dps = D + 40
    phi = (1+mp.sqrt(5))/2
    l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
    L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
    L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
    K5m = L5m + lf**5/120; K5q = L5q - mp.mpf(2)/15*lf**5
    K4m = L4m - lf**4/24;  K4q = L4q + lf**4/3
    K3m = L3m + lf**3/6
    atoms = {}
    for a in range(4):
        b = 3-a
        nm = 'z2' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
        atoms[nm] = z2*l2**a*lf**b
    for a in range(3):
        b = 2-a
        nm = 'z3' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
        atoms[nm] = z3*l2**a*lf**b
    atoms['z2^2l2'] = z2**2*l2; atoms['z2^2lf'] = z2**2*lf
    atoms['z5'] = z5; atoms['z2z3'] = z2*z3
    atoms['z2L3m'] = z2*L3m
    atoms['K5m'] = K5m; atoms['K5q'] = K5q
    atoms['K4ml2'] = K4m*l2; atoms['K4mlf'] = K4m*lf
    atoms['K4ql2'] = K4q*l2; atoms['K4qlf'] = K4q*lf
    atoms['K3ml2^2'] = K3m*l2**2; atoms['K3ml2lf'] = K3m*l2*lf; atoms['K3mlf^2'] = K3m*lf**2
    INV = ['z2l2^3','z2l2^2lf','z2l2lf^2','z2lf^3','z3l2^2','z3l2lf','z3lf^2','z2^2l2','z2^2lf',
           'z5','z2z3','z2L3m','K5m','K5q','K4ml2','K4mlf','K4ql2','K4qlf','K3ml2^2','K3ml2lf','K3mlf^2']
    kpool = pickle.load(open(os.path.join(HERE, 'kpool_1550.pkl'),'rb'))
    for nm, d in kpool.items():
        atoms[nm] = mp.mpf(d['K'])
        INV.append(nm)
    return atoms, INV

def main():
    D, TAB = int(sys.argv[1]), int(sys.argv[2])
    atoms, INV = setup(D)
    mp.mp.dps = D + 40
    r5 = mp.sqrt(5)
    vals = []
    for nm in INV:
        vals += [atoms[nm], r5*atoms[nm]]
    n = len(INV)
    C = mp.mpf(10)**int(D*0.92)
    rows = []
    for i in range(n):
        row = [0]*(n+2); row[i] = 1
        row[n] = int(mp.nint(C*vals[2*i])); row[n+1] = int(mp.nint(C*vals[2*i+1]))
        rows.append(row)
    M = fmpz_mat(rows).lll()
    print("LLL output shape:", M.nrows(), M.ncols())
    for i in range(min(4, n)):
        rv = [int(M[i,j]) for j in range(n+2)]
        t0 = sum(rv[j]*rows[j][n] for j in range(n)); t1 = sum(rv[j]*rows[j][n+1] for j in range(n+2) if j < n)
        print("row", i, "epart_nnz", sum(1 for x in rv[:n] if x), "tail", rv[n:], "expected", (t0, t1), "membership:", rv[n:][0]==t0 and rv[n:][1]==t1)
    rels = []
    for i in range(n):
        v = [int(M[i,j]) for j in range(n)]
        res1 = abs(sum(mp.mpf(v[j])*vals[2*j] for j in range(n)))
        res5 = abs(sum(mp.mpf(v[j])*vals[2*j+1] for j in range(n)))
        maxc = max(abs(x) for x in v)
        rels.append((max(mp.log10(res1+mp.mpf(10)**(-D-30)), mp.log10(res5+mp.mpf(10)**(-D-30))), maxc, v))
    rels.sort(key=lambda z: z[0])
    print("best relations (log10 residual, height):")
    for r0 in rels[:12]:
        print("  resid 1e%.1f  height %.3e" % (r0[0], r0[1]))
    # genuine relations: residual ~ 1e-0.9D
    gens = [v for (lg, mc, v) in rels if lg < -0.85*D]
    print("genuine relations found:", len(gens))
    # RREF over Q
    import sympy as sp
    if gens:
        Mk = sp.Matrix(gens)
        R, piv = Mk.rref()
        print("kernel rank:", len(piv), " pivot cols:", piv)
        for row in R.tolist():
            nz = {INV[j]: str(sp.Rational(c)) for j, c in enumerate(row) if c != 0}
            print("  REL:", nz)
    pickle.dump({'INV': INV, 'gens': gens}, open(os.path.join(HERE, 'kernel_gens.pkl'),'wb'))

if __name__ == "__main__":
    main()
