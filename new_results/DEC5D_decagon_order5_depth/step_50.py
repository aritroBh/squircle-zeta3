
gengine_dec = '''"""gengine_dec.py — numerical G(w;1) table for words over {0,1,-phi,phi^2}.
Same method as OCT5D2's validated gengine.py: exact-form power series at t=0 truncated at t^N,
evaluated at t=1/2; Chen concatenation to t=1. Caches to gtab_dec_d{dps}_N{N}.pkl."""
import os, math, pickle
import mpmath as mp

def build(dps=170, N=480):
    mp.mp.dps = dps
    phi = (1+mp.sqrt(5))/2
    L0, L1 = mp.mpf(0), mp.mpf(1)
    LP, LQ = -phi, phi**2
    LETTERS = [L0, L1, LP, LQ]
    inv = [mp.mpf(0)] + [mp.mpf(1)/j for j in range(1, N+2)]
    invp = {}
    for j in range(1, N+2):
        v = mp.mpf(1)/j; pw = v
        for k in range(1, 8):
            invp[(j,k)] = pw; pw *= v
    fact = [mp.mpf(math.factorial(i)) for i in range(8)]
    series = {(): [[mp.mpf(1)]] + [[mp.mpf(0)]]*N}
    words_by_len = [[()]]
    for length in range(1, 6):
        ws = []
        for w in words_by_len[length-1]:
            for a in LETTERS: ws.append((a,)+w)
        words_by_len.append(ws)
    for length in range(1, 6):
        L = length
        for w in words_by_len[length]:
            a, wr = w[0], w[1:]
            sw = series[wr]
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
    def mrev(v):
        return tuple(mp.mpf(1)-x for x in reversed(v))
    G1 = {}
    for w in Ghalf:
        tot = mp.mpf(0)
        for i in range(len(w)+1):
            u, v = w[:i], w[i:]
            tot += ((-1)**len(u)) * Ghalf[mrev(u)] * Ghalf[v]
        G1[w] = tot
    return G1, dict(L0=L0, L1=L1, LP=LP, LQ=LQ)

def gtab(dps=170, N=480, cachedir="."):
    fn = os.path.join(cachedir, f"gtab_dec_d{dps}_N{N}.pkl")
    if os.path.exists(fn):
        raw, letters = pickle.load(open(fn, "rb"))
        mp.mp.dps = dps
        return {k: mp.mpf(v) for k, v in raw.items()}, letters
    G1, letters = build(dps, N)
    raw = {k: mp.nstr(v, dps+10) for k, v in G1.items()}
    pickle.dump((raw, letters), open(fn, "wb"))
    return G1, letters
'''
open("gengine_dec.py","w").write(gengine_dec)
import gengine_dec
G1, letters = gengine_dec.gtab(120, 380, cachedir=".")
mp = __import__('mpmath').mp
mp.dps = 120
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
phi = (1+mp.sqrt(5))/2
# validations
r1 = G1[(L0,)*4+(LP,)] + mp.polylog(5, 1/LP)
print("G(0^4,-phi;1) + Li_5(-1/phi):", mp.nstr(abs(r1), 3))
r2 = G1[(L0,)*4+(LQ,)] + mp.polylog(5, 1/LQ)
print("G(0^4,phi^2;1) + Li_5(1/phi^2):", mp.nstr(abs(r2), 3))
gpq = mp.quad(lambda t: mp.log(1-t/LQ)/(t-LP), [0,1])
print("G(-phi,phi^2;1) vs quad:", mp.nstr(abs(G1[(LP,LQ)]-gpq), 3))
gqp = mp.quad(lambda t: mp.log(1-t/LP)/(t-LQ), [0,1])
print("G(phi^2,-phi;1) vs quad:", mp.nstr(abs(G1[(LQ,LP)]-gqp), 3))
sh = G1[(LP,)]*G1[(LQ,)] - G1[(LP,LQ)] - G1[(LQ,LP)]
print("shuffle G(p)G(q)=G(p,q)+G(q,p):", mp.nstr(abs(sh), 3))
r3 = G1[(L0,)*3+(LP,)] + mp.polylog(4, 1/LP)
print("G(0^3,-phi;1) + Li_4(-1/phi):", mp.nstr(abs(r3), 3))
r4 = G1[(LP,)] - mp.log(1-1/LP)
print("G(-phi;1) vs log(1+1/phi):", mp.nstr(abs(r4), 3))
