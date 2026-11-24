"""instr_deg.py (referee): reuse prove_adj78_fast's own functions to report, for (delta,s):
   - degrees (num-den) of every cancelled rational term, max over terms
   - degree of the lcm denominator DEN as the script builds it
   - bound on deg of the numerator polynomial N = (lhs - sum terms)*DEN :  maxdeg + degDEN
   - for each term, the polynomial P_j = num_j*(DEN/den_j)/content_j : degree and #terms (non-vacuity of the cancellation)
"""
import sys, time
import sympy as sp
import prove_adj78_fast as A
p = A.p
def run(delta, s, mod=A):
    t0 = time.time()
    lhs = sp.cancel(mod.eR(delta, s, mod.d, p) + mod.eR(delta, s+1, mod.d, p))
    terms = mod.feed_terms(delta, s)
    exprs = [lhs] + [-t for t in terms]
    rats = [sp.fraction(sp.cancel(e)) for e in exprs]
    fac_pow = {}; deninfo = []
    for num, den in rats:
        content, flist = sp.factor_list(den)
        cur = {}
        for f, e in flist:
            fm, lc = mod.monic(f)
            content *= lc**e
            if fm == 1: continue
            fac_pow[fm] = max(fac_pow.get(fm, 0), e)
            cur[fm] = e
        deninfo.append((content, cur))
    degs = [sp.degree(n, p) - sp.degree(dn, p) if n != 0 else None for n, dn in rats]
    degDEN = sum(sp.degree(f, p) * e for f, e in fac_pow.items())
    zero_terms = sum(1 for n, dn in rats if n == 0)
    maxdeg = max(x for x in degs if x is not None)
    sizes = []
    DEN = sp.Integer(1)
    for f, e in fac_pow.items(): DEN *= f**e
    DEN = sp.expand(DEN)
    for (num, den), (content, cur) in zip(rats, deninfo):
        quot = DEN
        for f, e in cur.items():
            for _ in range(e):
                q, r = sp.div(quot, f, p); assert r == 0; quot = q
        P = sp.Poly(sp.expand(num*quot/content), p)
        sizes.append((P.degree(), len(P.terms())))
    print(f'({delta},{s}): #exprs={len(rats)} (incl. lhs), identically-zero terms={zero_terms}, max term degree={maxdeg}, '
          f'#distinct lcm factors={len(fac_pow)}, deg DEN={degDEN}, => deg N <= {maxdeg+degDEN}; '
          f'per-term P_j (deg,#coeff) range: deg {min(x[0] for x in sizes)}..{max(x[0] for x in sizes)}, #coeff max {max(x[1] for x in sizes)}  ({round(time.time()-t0,1)}s)', flush=True)
    return maxdeg + degDEN
if __name__ == '__main__':
    run(int(sys.argv[1]), int(sys.argv[2]))
