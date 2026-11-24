"""assemble.py -- Task B: assemble c(w;alpha) = (1/(d-1)!) sum beta N from exact kernels + mpmath N-integrals."""
import pickle
from fractions import Fraction as Fr
from math import factorial
from mpmath import mp, mpf, quad, log

def qpart(t): return Fr(*t[0])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq

def load_kernels(name):
    import os
    A2 = '<project>/research/08_theory/ALPHA2_angle_leading_law'
    fn = f'{A2}/kern_{name}.pkl'
    if os.path.exists(f'{A2}/kernx_{name}.pkl'): fn = f'{A2}/kernx_{name}.pkl'
    data = pickle.load(open(fn,'rb'))
    return data

def beta_table(name, kmax=16):
    """extract beta_d^{(ij)}(rho) = leading coeff of residue polynomial, keys (d, (i,j), Fr pole)."""
    data = load_kernels(name)
    K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
    tab = {}
    for d in range(1, 8):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, kmax+1):
                P = K[k][key]
                for rho, cf in P['real'].items(): seqs.setdefault(rho, []).append(qpart(cf))
            for rho, s in seqs.items():
                if len(s) < d + 3: continue
                top = fdiff(list(s[1:]), d-1)
                assert all(x == top[0] for x in top)
                beta = top[0]/factorial(d-1)
                if beta != 0: tab[(d, key, Fr(*rho))] = beta
    return c, r, tab

def N_integrals(r, dmax=7, dps=55):
    """N_{(i,j),rho} = int_r^1 lt^i lmt^j dw/(w-rho), for i+j <= dmax, rho in {0,1,-1,-r}."""
    mp.dps = dps
    rm = mpf(r.numerator)/r.denominator
    def t(w): return (w-rm)*(w+1)/(2*(1-rm)*w)
    poles = [mpf(0), mpf(1), mpf(-1), -rm]
    N = {}
    for d in range(1, dmax+1):
        for i in range(d, (d-1)//2, -1):
            j = d - i
            for rho in poles:
                f = lambda w: log(t(w))**i * log(1-t(w))**j / (w - rho)
                N[((i,j), rho)] = quad(f, [rm, 1])
    return N

def assemble_cw(name, wmax=7, dps=55):
    """returns dict w -> c(w; alpha) as mpf, plus c, r."""
    c, r, tab = beta_table(name)
    N = N_integrals(r, dmax=wmax-1, dps=dps)
    rm = mpf(r.numerator)/r.denominator
    out = {}
    for w in range(2, wmax+1):
        d = w - 1
        tot = mpf(0)
        for (dd, key, rho), beta in tab.items():
            if dd != d: continue
            rhom = mpf(rho.numerator)/rho.denominator
            tot += mpf(beta.numerator)/beta.denominator * N[(key, rhom)]
        out[w] = tot/factorial(d-1)
    return out, c, r
