import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
import mpmath as mp
mp.mp.dps = 170
N = 480   # series truncation

# letters as mpf
r2 = mp.sqrt(2)
L0, L1 = mp.mpf(0), mp.mpf(1)
LP, LQ = 2+r2, -(1+r2)
LETTERS = [L0, L1, LP, LQ]

# ---- series at t=0, evaluated at t=1/2 ----
# series(w): list m=0..N of lists l=0..len(w) of mpf coeffs (t^m log^l t)
sys.setrecursionlimit(10000)
t0 = time.time()

# precompute 1/j^k tables
inv = [mp.mpf(0)] + [mp.mpf(1)/j for j in range(1, N+2)]
invp = {}  # invp[(j,k)] = 1/j^k
for j in range(1, N+2):
    v = mp.mpf(1)/j
    pw = v
    for k in range(1, 8):
        invp[(j,k)] = pw
        pw *= v

import math
fact = [mp.mpf(math.factorial(i)) for i in range(8)]

series = {(): [[mp.mpf(1)]]}
words_by_len = [[()]]
for length in range(1, 6):
    ws = []
    for w in words_by_len[length-1]:
        for a in LETTERS:
            ws.append((a,)+w)   # prepend: w = (a,)+w'  (w1 = a, suffix w')
    words_by_len.append(ws)

def integrate_series(sw, L):
    """int_0^t sw(s)/s ds ; sw[m][l] -> out with out[m][l], plus pure-log terms."""
    out = [[mp.mpf(0)]*(L+2) for _ in range(N+1)]
    for m in range(N+1):
        row = sw[m]
        if m == 0:
            for l in range(L+1):
                c = row[l]
                if c: out[0][l+1] += c*inv[l+1]
        else:
            for l in range(L+1):
                c = row[l]
                if not c: continue
                o = out[m]
                # t^m sum_{i<=l} (-1)^{l-i} l!/i! / m^{l-i+1} log^i
                sgn = mp.mpf(1)
                for i in range(l, -1, -1):
                    o[i] += c*sgn*fact[l]/fact[i]*invp[(m, l-i+1)]
                    sgn = -sgn
    return out

for length in range(1, 6):
    L = length
    for w in words_by_len[length]:
        a, wr = w[0], w[1:]
        sw = series[wr]
        if a == 0:
            out = integrate_series(sw, L-1)
        else:
            # multiply by -1/(s-a) via prefix recurrence M_m = c'_m + M_{m-1}/a ; product coeff = -M_m/a
            # then integrate termwise (no pure-log terms)
            out = [[mp.mpf(0)]*(L+1) for _ in range(N+1)]
            M = [mp.mpf(0)]*L      # M_{m-1} per log-degree (w' has degrees 0..L-1)
            for m in range(N+1):
                row = sw[m]
                for l in range(L):
                    c = row[l]
                    M[l] = (c + M[l])/a if (c or M[l]) else mp.mpf(0)
                if m == N: break
                om = out[m+1]
                j1 = m+1
                for l in range(L):
                    Mc = M[l]
                    if not Mc: continue
                    pc = -Mc/a
                    sgn = mp.mpf(1)
                    for i in range(l, -1, -1):
                        om[i] += pc*sgn*fact[l]/fact[i]*invp[(j1, l-i+1)]
                        sgn = -sgn
        series[w] = out
print("series built in %.1fs, words: %d" % (time.time()-t0, len(series)))

# evaluate at t=1/2
half = mp.mpf('0.5')
loghalf = -mp.log(2)
powh = [half**m for m in range(N+1)]
powl = [loghalf**l for l in range(7)]
t0 = time.time()
Ghalf = {}
for w, sw in series.items():
    tot = mp.mpf(0)
    for m in range(N+1):
        row = sw[m]
        pm = powh[m]
        for l in range(len(row)):
            c = row[l]
            if c: tot += c*pm*powl[l]
    Ghalf[w] = tot
print("evaluated in %.1fs" % (time.time()-t0))

# checks
print("G(0;1/2) - log(1/2) =", mp.nstr(Ghalf[(L0,)]-loghalf, 4))
print("G(1;1/2) - log(1/2) =", mp.nstr(Ghalf[(L1,)]-loghalf, 4))
print("G(p;1/2) - log(1-1/(2p)) =", mp.nstr(Ghalf[(LP,)]-mp.log(1-1/(2*LP)), 4))
print("G(0^5;1/2) - log^5(1/2)/120 =", mp.nstr(Ghalf[(L0,)*5]-loghalf**5/120, 4))
print("G(1^4;1/2) - log^4(1/2)/24 =", mp.nstr(Ghalf[(L1,)*4]-loghalf**4/24, 4))
print("G(0,q;1/2) vs quad:", mp.nstr(Ghalf[(L0,LQ)], 20), mp.nstr(mp.quad(lambda t: mp.log(1-t/LQ)/t, [0,half]), 20))
