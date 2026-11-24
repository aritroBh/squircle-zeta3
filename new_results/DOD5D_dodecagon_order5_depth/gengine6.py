"""gengine6.py -- numerical G(w;1) table for words over slots {0,1,2,3,4,5} =
{0, 1, p0, p1, p2, p3} (p_i = -1/kappa_i). Words are tuples of slot indices.
Same method as OCT5D2's gengine.py: exact-form power series at t=0 (truncated at t^N),
evaluated at t=1/2; Chen concatenation to t=1, mrev slot map 0<->1, 2<->3, 4<->5
(exact: 1-p0 = p1, 1-p1 = p0, 1-p2 = p3, 1-p3 = p2).
Lazy: builds series only for the closure of the needed words (suffixes + mrev prefixes).
Caches to gtab6_d{dps}_N{N}.pkl (values as mpf strings, keys slot tuples).
"""
import os, math, pickle
import mpmath as mp

MREV = (1, 0, 3, 2, 5, 4)
def letter_mpf():
    r3 = mp.sqrt(3)
    return [mp.mpf(0), mp.mpf(1), 4+2*r3, -(3+2*r3), -(1+r3)/2, (3+r3)/2]

def build(needed, dps=120, N=280):
    mp.mp.dps = dps
    LET = letter_mpf()
    words = set()
    for w in needed:
        for i in range(len(w)+1):
            words.add(w[i:])
        for i in range(len(w)+1):
            words.add(tuple(MREV[x] for x in reversed(w[:i])))
    words.add(())
    words_by_len = {}
    for w in words:
        words_by_len.setdefault(len(w), []).append(w)
    inv = [mp.mpf(0)] + [mp.mpf(1)/j for j in range(1, N+2)]
    invp = {}
    for j in range(1, N+2):
        v = mp.mpf(1)/j; pw = v
        for k in range(1, 8):
            invp[(j,k)] = pw; pw *= v
    fact = [mp.mpf(math.factorial(i)) for i in range(8)]
    series = {(): [[mp.mpf(1)]] + [[mp.mpf(0)]]*N}
    maxl = max(map(len, words))
    for length in range(1, maxl+1):
        for w in words_by_len.get(length, []):
            L = length
            a = LET[w[0]]
            sw = series[w[1:]]
            if a == 0:
                out = [[mp.mpf(0)]*(L+1) for _ in range(N+1)]
                for m in range(N+1):
                    row = sw[m]
                    if m == 0:
                        for l in range(L):
                            c = row[l]
                            if c: out[0][l+1] += c*inv[l+1]
                    else:
                        o = out[m]
                        for l in range(L):
                            c = row[l]
                            if not c: continue
                            sgn = mp.mpf(1)
                            for i in range(l, -1, -1):
                                o[i] += c*sgn*fact[l]/fact[i]*invp[(m, l-i+1)]
                                sgn = -sgn
            else:
                out = [[mp.mpf(0)]*(L+1) for _ in range(N+1)]
                M = [mp.mpf(0)]*L
                for m in range(N+1):
                    row = sw[m]
                    for l in range(L):
                        c = row[l]
                        if c or M[l]: M[l] = c + M[l]/a
                    if m == N: break
                    om = out[m+1]; j1 = m+1
                    for l in range(L):
                        Mc = M[l]
                        if not Mc: continue
                        pc = -Mc/a
                        sgn = mp.mpf(1)
                        for i in range(l, -1, -1):
                            om[i] += pc*sgn*fact[l]/fact[i]*invp[(j1, l-i+1)]
                            sgn = -sgn
            series[w] = out
    half = mp.mpf('0.5'); loghalf = -mp.log(2)
    powh = [half**m for m in range(N+1)]
    powl = [loghalf**l for l in range(7)]
    Ghalf = {}
    for w, sw in series.items():
        tot = mp.mpf(0)
        for m in range(N+1):
            row = sw[m]; pm = powh[m]
            for l in range(len(row)):
                c = row[l]
                if c: tot += c*pm*powl[l]
        Ghalf[w] = tot
    G1 = {(): mp.mpf(1)}
    for w in series:
        if w == (): continue
        tot = mp.mpf(0)
        for i in range(len(w)+1):
            u, v = w[:i], w[i:]
            mu = tuple(MREV[x] for x in reversed(u))
            tot += ((-1)**len(u)) * Ghalf[mu] * Ghalf[v]
        G1[w] = tot
    return G1

def gtab(needed, dps=120, N=280, cachedir="."):
    fn = os.path.join(cachedir, f"gtab6_d{dps}_N{N}.pkl")
    if os.path.exists(fn):
        raw = pickle.load(open(fn, "rb"))
        mp.mp.dps = dps
        return {k: mp.mpf(v) for k, v in raw.items()}, letter_mpf()
    G1 = build(needed, dps, N)
    raw = {k: mp.nstr(v, dps+10) for k, v in G1.items()}
    pickle.dump(raw, open(fn, "wb"))
    mp.mp.dps = dps
    return G1, letter_mpf()
