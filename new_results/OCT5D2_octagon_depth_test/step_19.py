import os, sys, time, pickle, math
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
import mpmath as mp
mp.mp.dps = 170
N = 480
r2 = mp.sqrt(2)
L0, L1 = mp.mpf(0), mp.mpf(1)
LP, LQ = 2+r2, -(1+r2)
LETTERS = [L0, L1, LP, LQ]

inv = [mp.mpf(0)] + [mp.mpf(1)/j for j in range(1, N+2)]
invp = {}
for j in range(1, N+2):
    v = mp.mpf(1)/j; pw = v
    for k in range(1, 8):
        invp[(j,k)] = pw; pw *= v
fact = [mp.mpf(math.factorial(i)) for i in range(8)]

def build_all(N, dps=None):
    series = {(): [[mp.mpf(1)]] + [[mp.mpf(0)]]*N}
    words_by_len = [[()]]
    for length in range(1, 6):
        ws = []
        for w in words_by_len[length-1]:
            for a in LETTERS: ws.append((a,)+w)
        words_by_len.append(ws)
    t0 = time.time()
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
                        if c or M[l]: M[l] = c + M[l]/a       # M_m = c'_m + M_{m-1}/a
                    if m == N: break
                    om = out[m+1]; j1 = m+1
                    for l in range(L):
                        Mc = M[l]
                        if not Mc: continue
                        pc = -Mc/a                            # product coeff at s^m
                        sgn = mp.mpf(1)
                        for i in range(l, -1, -1):
                            om[i] += pc*sgn*fact[l]/fact[i]*invp[(j1, l-i+1)]
                            sgn = -sgn
            series[w] = out
    return series, words_by_len

t0=time.time()
series, words_by_len = build_all(N)
print("series built in %.1fs" % (time.time()-t0))

half = mp.mpf('0.5'); loghalf = -mp.log(2)
powh = [half**m for m in range(N+1)]
powl = [loghalf**l for l in range(7)]
t0 = time.time()
Ghalf = {}
for w, sw in series.items():
    tot = mp.mpf(0)
    for m in range(N+1):
        row = sw[m]; pm = powh[m]
        for l in range(len(row)):
            c = row[l]
            if c: tot += c*pm*powl[l]
    Ghalf[w] = tot
print("evaluated in %.1fs" % (time.time()-t0))

print("G(p;1/2) - log(1-1/(2p)) =", mp.nstr(Ghalf[(LP,)]-mp.log(1-1/(2*LP)), 4))
print("G(q;1/2) - log(1-1/(2q)) =", mp.nstr(Ghalf[(LQ,)]-mp.log(1-1/(2*LQ)), 4))
print("G(0,q;1/2) + Li2(1/(2q)) =", mp.nstr(Ghalf[(L0,LQ)]+mp.polylog(2, 1/(2*LQ)), 4))
print("G(0,0,p;1/2) + Li3(1/(2p)) =", mp.nstr(Ghalf[(L0,L0,LP)]+mp.polylog(3, 1/(2*LP)), 4))
print("G(1;1/2)-log(1/2) =", mp.nstr(Ghalf[(L1,)]-loghalf, 4))
print("G(0^5;1/2)-log^5/120 =", mp.nstr(Ghalf[(L0,)*5]-loghalf**5/120, 4))
