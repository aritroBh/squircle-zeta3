import json, time
from mpmath import mp, mpf, log, pi, sqrt, exp, factorial, matrix, lu_solve, quad, nstr
t0=time.time()
mp.dps = 60
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
a_d = [mpf(x) for x in j120['d']]
b = [a_d[k]/factorial(k) for k in range(121)]   # Borel coefficients

def pade_borel_sum(b, d, L, M):
    # denominator: sum_{j=0}^M q_j b_{k-j} = 0, k = L+1..L+M, q_0=1
    A = matrix(M, M); Y = matrix(M, 1)
    for i in range(M):
        for jj in range(M):
            A[i, jj] = b[L - jj + i + 1 - 1] if L - jj + i + 1 - 1 >= 0 else mpf(0)
        Y[i] = -b[L + i + 1]
    qv = lu_solve(A, Y)
    q = [mpf(1)] + [qv[i] for i in range(M)]
    p = []
    for k in range(L + 1):
        s = mpf(0)
        for jj in range(min(k, M) + 1):
            s += q[jj] * b[k - jj]
        p.append(s)
    def P(z):
        num = sum(p[k] * z**k for k in range(L + 1))
        den = sum(q[jj] * z**jj for jj in range(M + 1))
        return num / den
    return quad(lambda t: exp(-t) * P(t * d), [0, mp.inf])

# my own direct quadrature reference (from step_23, recompute quickly at dps 45)
mp.dps = 45
def S_graph(p):
    p = mpf(p)
    def integrand(x, y):
        z = (1 - x**p - y**p)**(1/p)
        zx = -x**(p-1) * z**(1-p); zy = -y**(p-1) * z**(1-p)
        return sqrt(1 + zx**2 + zy**2)
    def inner(x):
        ymax = (1 - x**p)**(1/p)
        return ymax * quad(lambda v: integrand(x, ymax*v), [0, 1])
    return 8 * quad(inner, [0, 1])
from mpmath import log10
S1 = 4*sqrt(3)
for ds in ["0.01", "-0.01", "0.03", "0.1"]:
    d = mpf(ds)
    ref = S_graph(1 + d)/S1
    mp.dps = 60
    val = pade_borel_sum(b, d, 55, 55)
    dig = -log10(abs(val - ref)/abs(ref))
    print(f"d={ds}: Pade-Borel[55/55] vs own quadrature: {float(dig):.1f} digits  [{time.time()-t0:.0f}s]")
    mp.dps = 45
