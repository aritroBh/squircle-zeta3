
import mpmath as mp, json
mp.mp.dps = 80
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
d   = [mp.mpf(x) for x in h1['d']]
A3 = 2/mp.sqrt(mp.pi**2 + mp.log(2)**2)

def smoothed_root(a, lo, hi):
    """least squares fit of (|a_k|/k!)^(1/k) = A*(1+(b*log k + c)/k)"""
    ks = list(range(lo, hi+1))
    ys = [(abs(a[k])/mp.factorial(k))**(1/k) for k in ks]
    # fit y = A + A*b*log(k)/k + A*c/k  linear in (A, Ab, Ac)
    S = [[mp.mpf(0)]*3 for _ in range(3)]; T = [mp.mpf(0)]*3
    basis = lambda k: [mp.mpf(1), mp.log(k)/k, 1/k]
    for k, y in zip(ks, ys):
        b = basis(k)
        for i in range(3):
            T[i] += b[i]*y
            for j in range(3): S[i][j] += b[i]*b[j]
    M = mp.matrix(S); v = mp.matrix(T)
    x = mp.lu_solve(M, v)
    return x[0], x[1]/x[0], x[2]/x[0]

for lo in [30, 40, 50, 60]:
    At, bt, ct = smoothed_root(tau, lo, 79)
    Ad, bd, cd = smoothed_root(d, lo, 79)
    print(f"window {lo}..79: tau A={mp.nstr(At,9)} beta={mp.nstr(bt,6)} | d A={mp.nstr(Ad,9)} beta={mp.nstr(bd,6)}")
print("predicted A_tau =", mp.nstr(A3,9), " A_d =", mp.nstr(1/mp.mpf('0.7775'),5), "(=|1-t*|/|t*|, below)")

# exact mapped singularity in d: d* = t*/(1-t*)
ts = (-mp.log(2) + mp.pi*1j)/2
ds = ts/(1-ts)
print("t* =", mp.nstr(ts,12), " |t*| =", mp.nstr(abs(ts),12))
print("d* = t*/(1-t*) =", mp.nstr(ds,12), " |d*| =", mp.nstr(abs(ds),12), " 1/|d*| =", mp.nstr(1/abs(ds),9))
