"""RULQ script 0: the table J_{i,j} = int_0^oo s^j log^i(1+e^{-s}) ds  (i >= 1, j >= 0, i+j <= K-1).
Composite Gauss-Legendre on [0,S] with unit-2 panels (same route as PIPDIV pd1, which was checked there).
Usage: ru0_Jtable.py K DPS DEG S   ->  writes ru0_J_K{K}_dps{DPS}_deg{DEG}_S{S}.txt  (lines "i j value")
Check: J_{1,j} = j! (1-2^{-j-1}) zeta(j+2) is printed."""
import sys, time
import mpmath as mp
import mpmath.calculus.quadrature as Q

K = int(sys.argv[1]); mp.mp.dps = int(sys.argv[2]); DEG = int(sys.argv[3]); S = int(sys.argv[4])
t0 = time.time()
xs, ws = zip(*Q.GaussLegendre(mp.mp).calc_nodes(DEG, mp.mp.prec))
J = [[mp.mpf(0)] * K for _ in range(K + 1)]
for a0 in range(0, S, 2):
    for x, w in zip(xs, ws):
        s = a0 + 1 + x
        L = mp.log1p(mp.exp(-s))
        sp = [mp.mpf(1)]
        for j in range(1, K): sp.append(sp[-1] * s)
        Lr = w
        for i in range(1, K):
            Lr = Lr * L
            row = J[i]
            for j in range(0, K - i):
                row[j] += Lr * sp[j]
fn = 'ru0_J_K%d_dps%d_deg%d_S%d.txt' % (K, mp.mp.dps, DEG, S)
with open(fn, 'w') as f:
    for i in range(1, K):
        for j in range(0, K - i):
            f.write('%d %d %s\n' % (i, j, mp.nstr(J[i][j], mp.mp.dps)))
err = max(abs(J[1][j] / (mp.factorial(j) * (1 - mp.mpf(2) ** (-j - 1)) * mp.zeta(j + 2)) - 1) for j in range(K - 1))
print('max rel err of J_{1,j} vs closed form: %s' % mp.nstr(err, 3))
print('wrote %s in %.0fs' % (fn, time.time() - t0))
