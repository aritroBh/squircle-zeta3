
import mpmath as mp, json
mp.mp.dps = 90
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
d   = [mp.mpf(x) for x in h1['d']]
mine = json.load(open('SRF1_n3_K60.json'))
w = [mp.mpf(x) for x in mine['w']]

def borel_pade_poles(a, N):
    """[N/N] Pade of B(t)=sum (a_k/k!) t^k; return denominator roots"""
    b = [a[k]/mp.factorial(k) for k in range(2*N+1)]
    # solve for denominator q_1..q_N:  sum_{j=0}^{N} q_j b_{m-j} = 0, m=N+1..2N, q_0=1
    M = mp.matrix(N, N); v = mp.matrix(N, 1)
    for i in range(N):
        m = N+1+i
        for j in range(N):
            M[i,j] = b[m-1-j]
        v[i] = -b[m]
    q = mp.lu_solve(M, v)
    q = [mp.mpf(1)] + [q[j] for j in range(N)]
    # roots of Q(t) = sum q_j t^j
    coeffs = q[::-1]  # highest degree first
    try:
        roots = mp.polyroots(coeffs, maxsteps=200, error=False)
    except Exception:
        roots = []
    return roots

tedge = (-mp.log(2)+mp.pi*1j)/2
print("predicted t_edge =", mp.nstr(tedge,12), " |t|=", mp.nstr(abs(tedge),10))
for nm, a, N in [('tau', tau, 30), ('d', d, 30), ('w', w, 24)]:
    roots = borel_pade_poles(a, N)
    roots = sorted(roots, key=lambda z: abs(z))
    print(f"\n{nm}-series [{N}/{N}] Pade-Borel, nearest poles:")
    for r in roots[:4]:
        print("   ", mp.nstr(r, 10), " |r| =", mp.nstr(abs(r),8))
