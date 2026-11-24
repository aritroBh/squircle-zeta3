import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = '''"""gengine2.py — memory-efficient numerical G(w;1) table, words over {0,1,p,q}.

Same method as OCT5D2 gengine.py (series at t=1/2 + Chen concatenation), but stores the
full power series only for words of length <= KEEP (default 3); longer words are built one
at a time and their value at t=1/2 accumulated on the fly. Peak memory ~ O(words_len<=KEEP+1 * N).
Caches (Ghalf, G1) to gtab2_d{dps}_N{N}.pkl.
"""
import os, time, math, pickle
import mpmath as mp

KEEP = 3

def _letters():
    r2 = mp.sqrt(2)
    return [mp.mpf(0), mp.mpf(1), 2+r2, -(1+r2)]

def build(dps=170, N=560, keep=KEEP):
    mp.mp.dps = dps
    LETTERS = _letters()
    invp = {}
    for j in range(1, N+2):
        v = mp.mpf(1)/j; pw = v
        for k in range(1, 8):
            invp[(j,k)] = pw; pw *= v
    fact = [mp.mpf(math.factorial(i)) for i in range(8)]
    half = mp.mpf('0.5'); loghalf = -mp.log(2)
    powh = [half**m for m in range(N+1)]
    powl = [loghalf**l for l in range(7)]

    def series_extend(a, sw, L):
        """extend series of word wr (series sw, |wr|=L) to word (a,)+wr -> rows 0..N, len L+2."""
        if a == 0:
            out = [[mp.mpf(0)]*(L+2) for _ in range(N+1)]
            for m in range(N+1):
                row = sw[m]
                if m == 0:
                    for l in range(L+1):
                        c = row[l]
                        if c: out[0][l+1] += c*invp[(1, l+1)]
                else:
                    o = out[m]
                    for l in range(L+1):
                        c = row[l]
                        if not c: continue
                        sgn = mp.mpf(1)
                        for i in range(l, -1, -1):
                            o[i] += c*sgn*fact[l]/fact[i]*invp[(m, l-i+1)]
                            sgn = -sgn
            return out
        else:
            out = [[mp.mpf(0)]*(L+2) for _ in range(N+1)]
            M = [mp.mpf(0)]*(L+1)
            ainv = mp.mpf(1)/a
            for m in range(N+1):
                row = sw[m]
                for l in range(L+1):
                    c = row[l]
                    if c or M[l]: M[l] = c + M[l]*ainv
                if m == N: break
                om = out[m+1]; j1 = m+1
                for l in range(L+1):
                    Mc = M[l]
                    if not Mc: continue
                    pc = -Mc*ainv
                    sgn = mp.mpf(1)
                    for i in range(l, -1, -1):
                        om[i] += pc*sgn*fact[l]/fact[i]*invp[(j1, l-i+1)]
                        sgn = -sgn
            return out

    def value_at_half(sw):
        tot = mp.mpf(0)
        for m in range(N+1):
            row = sw[m]; pm = powh[m]
            for l in range(len(row)):
                c = row[l]
                if c: tot += c*pm*powl[l]
        return tot

    words_by_len = [[()]]
    for length in range(1, 6):
        ws = []
        for w in words_by_len[length-1]:
            for a in LETTERS: ws.append((a,)+w)
        words_by_len.append(ws)

    Ghalf = {(): mp.mpf(1)}
    # series kept: length <= keep all; length keep+1 all; longer: one at a time
    series_low = {(): [[mp.mpf(1)]] + [[mp.mpf(0)]]*N}   # series for words of length <= keep
    parent = series_low          # series dict for words of current parent length
    for length in range(1, 6):
        cur = {}
        store = (length <= keep+1)
        for w in words_by_len[length]:
            a, wr = w[0], w[1:]
            sw = series_extend(a, parent[wr], length-1)
            if length <= keep+1:
                Ghalf[w] = value_at_half(sw)
                cur[w] = sw
            else:
                Ghalf[w] = value_at_half(sw)
                del sw
        if length <= keep:
            series_low.update(cur)
        if length == keep+1:
            parent = cur           # length keep+1 series, needed for keep+2 children
        elif length > keep+1:
            parent = None
        if length <= keep:
            parent = series_low
    # Chen
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
    # self-check: G(0^4,p;1) = -Li5(1/p)
    L0, LP = letters['L0'], letters['LP']
    err = abs(G1[(L0,)*4+(LP,)] + mp.polylog(5, 1/LP))
    print("G(0^4,p;1)+Li5(1/p) residual:", mp.nstr(err, 3))
'''
open("gengine2.py","w").write(src)
print("written gengine2.py")
