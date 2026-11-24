"""Write genfun.json (exact Newton data for every key/pole, d<=8) and genfun_d1_d3.md (readable tables for d=1,2,3).
G_r(eps) = c_r(d) eps^d + sum_{m=1}^d D_m eps^{d+m}/(1-eps)^m = eps^d [c_r(d) + sum D_m x^m],  x = eps/(1-eps),
c_r(k) = sum_{m=1}^d D_m binom(k-d-1, m-1) for k>=d+1 (D_m = Delta^{m-1} c_r at k=d+1)."""
import sys, os, pickle, glob, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
import sympy as sp
eps, x = sp.symbols('eps x')
S2 = sp.sqrt(2)
kern = {int(os.path.basename(f)[5:-4]): pickle.load(open(f, 'rb')) for f in glob.glob(os.path.join(HERE, "kern_*.pkl"))}
ks = sorted(kern); K = ks[-1]
POLE = {(0, 0): "0", (1, 1): "a", (-1, -1): "-a", (-1, 1): "b", (1, -1): "-b"}   # a=1+sqrt2, b=sqrt2-1 (so -b = 1-sqrt2)
def c(k, key, r):
    P = kern.get(k, {}).get(key)
    if P is None or r not in P['real']: return (Fr(0), Fr(0))
    z = P['real'][r]; return (z[0], z[1])
def sp_of(t): return sp.Rational(t[0].numerator, t[0].denominator) + sp.Rational(t[1].numerator, t[1].denominator)*S2
def fd(seq, m):
    s = list(seq)
    for _ in range(m): s = [(a2[0]-a1[0], a2[1]-a1[1]) for a1, a2 in zip(s, s[1:])]
    return s
out = {}; lines = []
DMAX = 8
for d in range(1, DMAX+1):
    keys = sorted({key for k in ks for key in kern[k] if sum(key) == d}, key=lambda kk: -kk[0])
    for key in keys:
        poles = sorted({r for k in ks for r in kern[k].get(key, {'real': {}})['real']})
        rec = {}
        for r in poles:
            seq = [c(k, key, r) for k in range(d+1, K+1)]
            D = []; s = list(seq)
            for m in range(1, d+1):
                D.append(s[0]); s = [(a2[0]-a1[0], a2[1]-a1[1]) for a1, a2 in zip(s, s[1:])]
            assert all(t == (0, 0) for t in s), ("NOT polynomial of degree <= d-1", key, r)
            rec[POLE[r]] = {"c_d": [str(c(d, key, r)[0]), str(c(d, key, r)[1])], "D": [[str(t[0]), str(t[1])] for t in D],
                            "n_independent_checks": len(seq) - d}
        out[f"{d}|{key[0]},{key[1]}"] = rec
        if d <= 3:
            lines.append(f"\n**key ({key[0]},{key[1]}), d={d}**  (poles v = 0, a, -a, b, -b with a=1+sqrt2, b=sqrt2-1; no pole at +b ever occurs)\n")
            lines.append("| pole | G_r(eps) = numerator/(1-eps)^m (expanded; sqrt2 explicit) | eps^d[c_d + D_1 x + ... + D_d x^d] data (c_d; D_1..D_d) |")
            lines.append("|---|---|---|")
            for r in poles:
                cd = sp_of(c(d, key, r))
                Ds = [sp_of(t) for t in [t_ for t_ in [ (rec[POLE[r]]['D'][m][0], rec[POLE[r]]['D'][m][1]) for m in range(d)]]] if False else None
                Dl = [sp.Rational(Fr(rec[POLE[r]]['D'][m][0]).numerator, Fr(rec[POLE[r]]['D'][m][0]).denominator) + sp.Rational(Fr(rec[POLE[r]]['D'][m][1]).numerator, Fr(rec[POLE[r]]['D'][m][1]).denominator)*S2 for m in range(d)]
                G = cd*eps**d + sum(Dl[m-1]*eps**(d+m)/(1-eps)**m for m in range(1, d+1))
                Gc = sp.cancel(sp.together(G))
                num_, den_ = sp.fraction(Gc)
                num_ = sp.expand(num_); den_ = sp.expand(den_)
                # normalise the denominator to (1-eps)^m times a rational number
                m_ = sp.degree(den_, eps)
                lead = sp.Poly(den_, eps).LC()
                sgn = lead * (-1)**m_          # den_ = lead*eps^m+... = sgn*(1-eps)^m
                num_n = sp.expand(num_ / sgn)
                lines.append(f"| {POLE[r]} | [{sp.sstr(sp.collect(num_n, eps))}] / (1-eps)^{m_} | c_d={sp.sstr(sp.expand(cd))}; D={[sp.sstr(sp.expand(t)) for t in Dl]} |")
                chk = sp.simplify(num_n/(1-eps)**m_ - G)
                assert chk == 0, ("normalisation error", key, r)
json.dump(out, open(os.path.join(HERE, "genfun.json"), "w"), indent=1)
open(os.path.join(HERE, "genfun_d1_d3.md"), "w").write("\n".join(lines) + "\n")
print("wrote genfun.json (", len(out), "keys, d<=8 ) and genfun_d1_d3.md; every sequence checked: Delta^d c_r(k) = 0 for k>=d+1 up to k =", K)
