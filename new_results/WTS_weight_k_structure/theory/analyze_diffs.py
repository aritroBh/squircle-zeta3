"""TASK 1: exact test of the conjecture on the kernels kern_k.pkl (all keys, all k available).
For every key (i,j), d=i+j, every real pole r, the residue c_r(k) in Q(sqrt2) is a sequence in k.
Test: forward differences of order d (in k, starting at k=d+1) vanish identically, i.e. c_r(k) is a
polynomial in k of degree <= d-1 for k >= d+1.  Report also the minimal polynomial degree seen, and
whether k=d belongs to the polynomial (it should not: exceptional).
Everything is exact (Fractions in Q(sqrt2): pairs).  Prints PASS/FAIL per (key).
"""
import sys, os, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
from math import comb

kern = {}
for fn in sorted(glob.glob(os.path.join(HERE, "kern_*.pkl"))):
    k = int(os.path.basename(fn)[5:-4])
    kern[k] = pickle.load(open(fn, 'rb'))
ks = sorted(kern)
print("levels available:", ks)

def realpart(P):
    return P

# residue table: kern[k][key]['real'] = {(r0,r1): (c0,c1,c2,c3)}  (r = r0 + r1 sqrt2, c likewise)
def poles_of(key):
    S = set()
    for k in ks:
        if key in kern[k]:
            S |= set(kern[k][key]['real'].keys())
    return sorted(S)

def coeff(k, key, r):
    P = kern[k].get(key)
    if P is None: return (Fr(0),)*4
    return P['real'].get(r, (Fr(0),)*4)

# check no +-i poles, no constants, c real
for k in ks:
    for key, P in kern[k].items():
        assert not any(any(x != 0 for x in e) for e in P['ipos']), (k, key, "ipos")
        assert not any(any(x != 0 for x in e) for e in P['ineg']), (k, key, "ineg")
        assert all(x == 0 for x in P['const']), (k, key, "const")
        for r, c in P['real'].items():
            assert c[2] == 0 and c[3] == 0, (k, key, r, c, "nonreal residue")
print("all kernels: no +-i poles, no constant part, real residues: OK")

def fdiff(seq, order):
    s = list(seq)
    for _ in range(order):
        s = [s[i+1] - s[i] for i in range(len(s)-1)]
    return s

allkeys = sorted({key for k in ks for key in kern[k]}, key=lambda kk: (sum(kk), -kk[0]))
summary = {}
for key in allkeys:
    d = sum(key)
    krange = [k for k in ks if k >= d+1]
    if len(krange) < d+1:
        summary[key] = ("NOTESTABLE", len(krange))
        continue
    poles = poles_of(key)
    # for each (pole, component 0/1) sequence
    okall = True; maxdeg = -1; ncomp = 0
    for r in poles:
        for comp in (0, 1):
            seq = [coeff(k, key, r)[comp] for k in krange]
            # minimal n such that order-(n+1) difference vanishes
            n = None
            for m in range(0, len(seq)):
                if all(x == 0 for x in fdiff(seq, m)):
                    n = m - 1; break
            ncomp += 1
            if n is None or n > d-1:
                okall = False
            if n is not None: maxdeg = max(maxdeg, n)
    nchecks = len(krange) - d   # number of independent d-th differences tested per sequence
    summary[key] = ("PASS" if okall else "FAIL", maxdeg, d-1, len(krange), nchecks)
    # does the polynomial extend to k=d ?
    ext = None
    if d in ks and key in kern[d]:
        ext = True
        for r in poles:
            for comp in (0, 1):
                seq = [coeff(k, key, r)[comp] for k in krange]
                # Lagrange extrapolate polynomial through first d points to k=d
                pts = list(zip(krange[:d], seq[:d]))
                val = Fr(0)
                for a, (ka, ya) in enumerate(pts):
                    term = ya
                    for b, (kb, yb) in enumerate(pts):
                        if a != b: term *= Fr(d - kb, ka - kb)
                    val += term
                if val != coeff(d, key, r)[comp]: ext = False
    summary[key] = summary[key] + (("k=d on polynomial: %s" % ext),)
for key in allkeys:
    print(key, summary[key])
