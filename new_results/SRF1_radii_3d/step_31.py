
import mpmath as mp, json
mp.mp.dps = 160
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
mine = json.load(open('SRF1_n3_K60.json'))
w = [mp.mpf(x) for x in mine['w']]
A3 = 2/mp.sqrt(mp.pi**2 + mp.log(2)**2)
tedge = (-mp.log(2)+mp.pi*1j)/2

def borel_pade_poles(a, N):
    b = [a[k]/mp.factorial(k)*A3**(-k) for k in range(2*N+1)]  # scaled: poles near |z|=1
    M = mp.matrix(N, N); v = mp.matrix(N, 1)
    for i in range(N):
        m = N+1+i
        for j in range(N): M[i,j] = b[m-1-j]
        v[i] = -b[m]
    q = mp.lu_solve(M, v)
    q = [mp.mpf(1)] + [q[j] for j in range(N)]
    roots = mp.polyroots(q[::-1], maxsteps=400, error=False)
    return [z/A3 for z in roots]  # unscale: t = z/A3... wait

# careful: if B(t) = sum b_k (t*A3)^k... define z = t*A3, series in z has poles at z* = t*A3 (|z*|=1)
def borel_pade_poles2(a, N):
    b = [a[k]/mp.factorial(k)*A3**(-k) for k in range(2*N+1)]
    M = mp.matrix(N, N); v = mp.matrix(N, 1)
    for i in range(N):
        m = N+1+i
        for j in range(N): M[i,j] = b[m-1-j]
        v[i] = -b[m]
    q = mp.lu_solve(M, v); q = [mp.mpf(1)]+[q[j] for j in range(N)]
    roots = mp.polyroots(q[::-1], maxsteps=400, error=False)
    return [z/A3 for z in roots]  # t = z/A3

print("predicted t_edge =", mp.nstr(tedge,12))
for nm, a, N in [('tau', tau, 30), ('w', w, 22)]:
    roots = borel_pade_poles2(a, N)
    roots = sorted(roots, key=lambda z: abs(z))
    print(f"\n{nm} [{N}/{N}] nearest Pade-Borel poles (t-plane):")
    for r in roots[:6]:
        print("   ", mp.nstr(r, 9), "  |r| =", mp.nstr(abs(r),7))
