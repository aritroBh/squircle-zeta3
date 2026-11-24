"""gengine2.py — memory-efficient numerical G(w;1) table, words over {0,1,p,q}.

Method of OCT5D2 gengine.py (exact-form series at t=1/2 + Chen concatenation), rewritten to
store full series only for words of length <= keep+1; longer words are built one at a time
and evaluated at t=1/2 on the fly. Optional multiprocessing over top-level words.
Caches (G1) to gtab2_d{dps}_N{N}.pkl (values as mpf strings).
"""
import os, time, math, pickle
import mpmath as mp

KEEP = 3

def _setup(dps, N):
    mp.mp.dps = dps
    r2 = mp.sqrt(2)
    LETTERS = [mp.mpf(0), mp.mpf(1), 2+r2, -(1+r2)]
    INVP = []                      # INVP[j][k] = (1/j)^k
    INVP.append([mp.mpf(0)]*9)
    for j in range(1, N+2):
        v = mp.mpf(1)/j; pw = v; row = [mp.mpf(0)]
        for k in range(1, 9):
            row.append(pw); pw *= v
        INVP.append(row)
    FRAT = [[mp.mpf(math.factorial(l))/math.factorial(i) for i in range(8)] for l in range(8)]
    half = mp.mpf('0.5'); loghalf = -mp.log(2)
    POWH = [half**m for m in range(N+1)]
    POWL = [loghalf**l for l in range(8)]
    return LETTERS, INVP, FRAT, POWH, POWL

def _extend(a, sw, L, N, INVP, FRAT):
    """series of (a,)+wr from series sw of wr (|wr|=L) -> rows 0..N, length L+2."""
    zero = mp.mpf(0)
    out = [[zero]*(L+2) for _ in range(N+1)]
    if a == 0:
        row0 = sw[0]; o0 = out[0]
        for l in range(L+1):
            c = row0[l]
            if c: o0[l+1] = c/(l+1)
        for m in range(1, N+1):
            row = sw[m]; o = out[m]; ip = INVP[m]
            for l in range(L+1):
                c = row[l]
                if not c: continue
                sgn = mp.mpf(1); fr = FRAT[l]
                for i in range(l, -1, -1):
                    o[i] += c*sgn*fr[i]*ip[l-i+1]
                    sgn = -sgn
        return out
    else:
        ainv = mp.mpf(1)/a
        M = [mp.mpf(0)]*(L+1)
        for m in range(N+1):
            row = sw[m]
            for l in range(L+1):
                c = row[l]
                if c or M[l]: M[l] = c + M[l]*ainv
            if m == N: break
            om = out[m+1]; ip = INVP[m+1]
            for l in range(L+1):
                Mc = M[l]
                if not Mc: continue
                pc = -Mc*ainv; sgn = mp.mpf(1); fr = FRAT[l]
                for i in range(l, -1, -1):
                    om[i] += pc*sgn*fr[i]*ip[l-i+1]
                    sgn = -sgn
        return out

def _value_at_half(sw, POWH, POWL):
    tot = mp.mpf(0)
    for m in range(len(sw)-1, -1, -1):
        row = sw[m]; pm = POWH[m]
        for l in range(len(row)):
            c = row[l]
            if c: tot += c*pm*POWL[l]
    return tot

def _words(LETTERS, maxlen=5):
    words_by_len = [[()]]
    for length in range(1, maxlen+1):
        ws = []
        for w in words_by_len[length-1]:
            for a in LETTERS: ws.append((a,)+w)
        words_by_len.append(ws)
    return words_by_len

def build(dps=170, N=560, keep=KEEP):
    mp.mp.dps = dps
    LETTERS, INVP, FRAT, POWH, POWL = _setup(dps, N)
    words_by_len = _words(LETTERS)
    Ghalf = {(): mp.mpf(1)}
    series_low = {(): [[mp.mpf(1)]] + [[mp.mpf(0)]]*N}
    parent = series_low
    for length in range(1, 6):
        cur = {}
        if length <= keep+1:
            for w in words_by_len[length]:
                sw = _extend(w[0], parent[w[1:]], length-1, N, INVP, FRAT)
                Ghalf[w] = _value_at_half(sw, POWH, POWL)
                cur[w] = sw
            if length <= keep:
                series_low.update(cur); parent = series_low
            else:
                parent = cur
        else:
            for w in words_by_len[length]:
                sw = _extend(w[0], parent[w[1:]], length-1, N, INVP, FRAT)
                Ghalf[w] = _value_at_half(sw, POWH, POWL)
                del sw
    def mrev(v):
        return tuple(mp.mpf(1)-x for x in reversed(v))
    G1 = {}
    for w, gw in Ghalf.items():
        tot = mp.mpf(0)
        for i in range(len(w)+1):
            u, v = w[:i], w[i:]
            tot += ((-1)**len(u)) * Ghalf[mrev(u)] * Ghalf[v]
        G1[w] = tot
    letters = dict(L0=LETTERS[0], L1=LETTERS[1], LP=LETTERS[2], LQ=LETTERS[3])
    return G1, letters

def gtab(dps=170, N=560, cachedir="."):
    fn = os.path.join(cachedir, f"gtab2_d{dps}_N{N}.pkl")
    if os.path.exists(fn):
        raw, letters = pickle.load(open(fn, "rb"))
        mp.mp.dps = dps
        return {k: mp.mpf(v) for k, v in raw.items()}, letters
    G1, letters = build(dps, N)
    raw = {k: mp.nstr(v, dps+10) for k, v in G1.items()}
    pickle.dump((raw, letters), open(fn, "wb"))
    return G1, letters

if __name__ == "__main__":
    import sys
    dps, N = int(sys.argv[1]), int(sys.argv[2])
    t0 = time.time()
    G1, letters = gtab(dps, N)
    print("built dps=%d N=%d in %.1fs, %d words" % (dps, N, time.time()-t0, len(G1)))
    L0, LP = letters['L0'], letters['LP']
    err = abs(G1[(L0,)*4+(LP,)] + mp.polylog(5, 1/LP))
    print("G(0^4,p;1)+Li5(1/p) residual:", mp.nstr(err, 3))
