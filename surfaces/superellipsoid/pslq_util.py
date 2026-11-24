# pslq_util.py -- integer-relation helper (mpmath.pslq) with residual reporting.
import mpmath as mp

def find_relation(target, basis, names, dps, maxcoeff=10**8, maxsteps=10**6, tol_digits=None):
    """Return (coeffs, residual) for an integer relation c0*target + sum ci*basis_i = 0, or None."""
    with mp.workdps(dps):
        vec = [mp.mpf(target)] + [mp.mpf(b) for b in basis]
        tol = mp.mpf(10)**(-(tol_digits if tol_digits else int(dps*0.9)))
        rel = mp.pslq(vec, tol=tol, maxcoeff=maxcoeff, maxsteps=maxsteps)
        if rel is None:
            return None
        res = mp.fsum(c*v for c, v in zip(rel, vec))
        return rel, res

def show(rel, names):
    c0 = rel[0]
    terms = []
    for c, n in zip(rel[1:], names):
        if c != 0:
            terms.append("%s*%s" % (mp.nstr(mp.mpf(-c)/c0, 12) if False else "(%d/%d)" % (-c, c0), n))
    return " + ".join(terms)
