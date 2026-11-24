"""RULQR script 2: Taylor coefficients of G(h) = Ehat(h)/h^3 from the Chebyshev-node values of rr1_values.py,
then N_k = (-log2)^{k+1}/(k+1)! + [h^{k-2}]G.   Usage: rr2_extract.py VALFILE H KMAX OUT [SPOW]"""
import sys, mpmath as mp
vf, H, KMAX, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
SPOW = int(sys.argv[5]) if len(sys.argv) > 5 else 3   # G = Ehat/h^SPOW
mp.mp.dps = 700
H = mp.mpf(H)
vals = []
for line in open(vf):
    g = line.split(' ', 1)[1].strip().strip('[]').split('+/-')[0]
    vals.append(mp.mpf(g))
D = len(vals)
# Chebyshev coefficients (nodes x_j = cos(pi(j+1/2)/D) on [-1,1], h = H(1+x)/2)
c = [2 * mp.fsum(vals[j] * mp.cos(m * mp.pi * (j + mp.mpf(1) / 2) / D) for j in range(D)) / D for m in range(D)]
c[0] /= 2
# monomial coefficients in h: T_m(2h/H - 1) via recurrence on coefficient lists
x = [mp.mpf(-1), 2 / H]
Tm1, T = [mp.mpf(1)], x[:]
G = [mp.mpf(0)] * (D + 1)
G[0] += c[0]
for k in range(2): G[k] += c[1] * T[k]
for m in range(2, D):
    Tn = [mp.mpf(0)] * (m + 1)
    for k, a in enumerate(T):
        Tn[k] += -2 * a; Tn[k + 1] += 2 * a * (2 / H)
    for k, a in enumerate(Tm1): Tn[k] -= a
    Tm1, T = T, Tn
    for k in range(m + 1): G[k] += c[m] * T[k]
L2 = mp.log(2)
with open(out, 'w') as fo:
    for k in range(KMAX + 1):
        N = (-L2) ** (k + 1) / mp.factorial(k + 1) + (G[k + 1 - SPOW] if k + 1 >= SPOW else 0)
        fo.write(mp.nstr(N, 60) + '\n')
print('tail |c_m| (last 3):', [mp.nstr(abs(t), 3) for t in c[-3:]])
