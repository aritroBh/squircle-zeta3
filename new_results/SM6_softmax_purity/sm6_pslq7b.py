import mpmath as mp, pickle, itertools
mp.mp.dps = 420
z = mp.zeta; L2 = mp.log(2)
CAND, cv, gams = pickle.load(open("sm6_w7_raw.pkl", "rb"))
g33, g15 = gams
Li = lambda n: mp.polylog(n, mp.mpf(1)/2)
Z6 = mp.quad(lambda x: (-mp.log(x))**4*mp.log(1-x)/(1+x), [0, mp.mpf(1)/2, 1], maxdegree=8)/24
# generators and weights
def build(X):
    gens = [("L", 1, L2), ("z2", 2, z(2)), ("z3", 3, z(3)), ("Li4", 4, Li(4)), ("z5", 5, z(5)), ("Li5", 5, Li(5)),
            ("Li6", 6, Li(6)), ("Z6", 6, Z6), ("z7", 7, z(7)), ("Li7", 7, Li(7)), ("X1", 7, X[0]), ("X2", 7, X[1])]
    mons = []
    def rec(i, w, name, val):
        if w == 7: mons.append((name, val)); return
        if i == len(gens): return
        n, wt, v = gens[i]
        k = 0
        while w + k*wt <= 7:
            rec(i+1, w+k*wt, name + ("*%s^%d" % (n, k) if k else ""), val*v**k)
            k += 1
    rec(0, 0, "", mp.mpf(1))
    return mons
if __name__ == "__main__":
    best = None
    for i, j in itertools.combinations(range(len(CAND)+1), 2):
        X = (cv[i], cv[j])
        mons = build(X)
        assert len(mons) == 21, len(mons)
        mp.mp.dps = 380
        r = mp.pslq([v for _, v in mons], maxcoeff=10**12, maxsteps=10**6)
        mp.mp.dps = 420
        print("pair", i, j, "independent" if r is None else "DEPENDENT", flush=True)
        if r is None:
            rr = mp.pslq([g33] + [v for _, v in mons], maxcoeff=10**14, maxsteps=10**6) if False else None
            mp.mp.dps = 380
            rr = mp.pslq([g33] + [v for _, v in mons], maxcoeff=10**14, maxsteps=10**6)
            mp.mp.dps = 420
            print("   gamma33 relation:", rr, flush=True)
            if rr: best = (i, j, rr, mons); break
    if best:
        i, j, rr, mons = best
        print("candidates used:", CAND[i] if i < len(CAND) else "(4,1,0,1+x) -> int (-log x)^4 log(1-x)/(1+x)", CAND[j] if j < len(CAND) else "extra")
        print("coefficient of gamma33:", rr[0])
        for c, (n, _) in zip(rr[1:], mons):
            if c: print("  %d  %s" % (c, n))
        pickle.dump((i, j, rr, [n for n, _ in mons]), open("sm6_w7_rel.pkl", "wb"))
