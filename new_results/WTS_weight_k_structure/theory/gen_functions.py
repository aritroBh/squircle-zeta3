"""TASK 1b: exact eps-generating functions  G_key(eps; v) = sum_k eps^k K^{(k)}_key(v)  for each key (i,j), d=i+j<=DMAX.
K^{(k)}(v) = sum_r c_r(k)/(v-r);  conjecture: c_r(k) polynomial in k of degree <= d-1 for k>=d+1.
Newton form: c_r(k) = sum_{m=1}^{d} D_m binom(k-d-1, m-1),  D_m = Delta^{m-1} c_r (k=d+1)   (k>=d+1)
=> G_r(eps) = c_r(d) eps^d + sum_m D_m eps^{d+m}/(1-eps)^m = eps^d [ c_r(d) + sum_m D_m x^m ],  x = eps/(1-eps).
The script RE-EXPANDS G_r(eps) as a series and compares with the data for every available k (2..15), exactly.
"""
import sys, os, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import sympy as sp
from fractions import Fraction as Fr
eps, x, vv = sp.symbols('eps x v')
S2 = sp.sqrt(2)
kern = {}
for fn in sorted(glob.glob(os.path.join(HERE, "kern_*.pkl"))):
    kern[int(os.path.basename(fn)[5:-4])] = pickle.load(open(fn, 'rb'))
ks = sorted(kern)
def q2sp(p0, p1):
    return sp.Rational(p0.numerator, p0.denominator) + sp.Rational(p1.numerator, p1.denominator)*S2
def pole_sp(r):
    return q2sp(r[0], r[1])
def coeff(k, key, r):
    P = kern[k].get(key)
    if P is None or r not in P['real']: return (Fr(0), Fr(0))
    c = P['real'][r]; return (c[0], c[1])
def cspx(k, key, r):
    c = coeff(k, key, r); return q2sp(*c)
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
out_lines = []
ALL_OK = True
for d in range(1, DMAX+1):
    keys = sorted({key for k in ks for key in kern[k] if sum(key) == d}, key=lambda kk: -kk[0])
    for key in keys:
        poles = sorted({r for k in ks for r in kern[k].get(key, {'real': {}})['real']})
        Gtot = {}
        for r in poles:
            seq = [cspx(k, key, r) for k in range(d+1, ks[-1]+1)]
            # Newton forward differences at k=d+1
            D = []; s = list(seq)
            for m in range(1, d+1):
                D.append(sp.expand(s[0]))
                s = [sp.expand(s[i+1]-s[i]) for i in range(len(s)-1)]
            # after d differences the remaining must be zero
            assert all(sp.simplify(t) == 0 for t in s), ("difference not zero", key, r)
            c_d = cspx(d, key, r)
            Gr = sp.expand(c_d*eps**d) + sum(D[m-1]*eps**(d+m)/(1-eps)**m for m in range(1, d+1))
            Gr = sp.together(Gr)
            num, den = sp.fraction(sp.cancel(Gr))
            num = sp.expand(num); 
            # verify series against data for all k
            ser = sp.series(Gr, eps, 0, ks[-1]+1).removeO()
            for k in ks:
                lhs = sp.expand(ser.coeff(eps, k)); rhs = sp.expand(cspx(k, key, r))
                if sp.simplify(lhs - rhs) != 0:
                    ALL_OK = False
                    print("MISMATCH", key, r, k, lhs, rhs)
            # and zero below level d
            for k in range(0, d):
                assert sp.simplify(ser.coeff(eps, k)) == 0
            Gtot[r] = (sp.factor(num, extension=[S2]), sp.factor(den), D, c_d)
        out_lines.append((key, poles, Gtot))
        print(f"=== key {key} (d={d}) ===")
        for r in poles:
            num, den, D, c_d = Gtot[r]
            print(f"  pole v = {sp.expand(pole_sp(r))}:  G = ({sp.expand(num)}) / ({den})")
            print(f"      c(d)={sp.expand(c_d)};  Newton D_m = {[sp.expand(t) for t in D]}")
print("ALL GENERATING FUNCTIONS REPRODUCE THE DATA EXACTLY (all k=2..%d, incl. k=d, zero below):" % ks[-1], ALL_OK)
pickle.dump({str(k): (k_, [ (r, [str(t) for t in Gtot[r][2]], str(Gtot[r][3]), str(Gtot[r][0]), str(Gtot[r][1])) for r in Gtot]) for (k, (k_, _, Gtot)) in enumerate(out_lines)}, open(os.path.join(HERE, 'genfun_summary.pkl'), 'wb'))
