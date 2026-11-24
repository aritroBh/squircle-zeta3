
"""
SRF3 step_05: Borel summability of the octahedron-end series along the real axis (Task 4).
The Borel singularities z0=(-log2+-i pi)/2 have Re z0 < 0, so the Laplace integral along R+
avoids them: the series of S(1+d)/S(1) in tau (and d) should be Borel-summable.
Test: Padé-Borel resummation from K=120 coefficients vs H1's 40-digit direct quadrature.
"""
import mpmath as mp, json, os
mp.mp.dps = 80
wd = os.path.dirname(os.path.abspath(__file__))
base = "<project>"
h = json.load(open(wd+"/SRF3_hi_coeffs.json"))
tau = [mp.mpf(x) for x in h['eta']]   # a_k = [tau^k] S/S(1)
dco = [mp.mpf(x) for x in h['d']]
K = h['K']
direct = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/direct_S_dps40.json")))
S1 = 4*mp.sqrt(3)

def pade(c, L, M):
    A = mp.matrix(M,M); r = mp.matrix(M,1)
    for i in range(M):
        for j in range(M):
            A[i,j] = c[L-M+i+j+1]
        r[i] = -c[L+i+1]
    q = mp.lu_solve(A,r); q = [mp.mpf(1)]+[q[i] for i in range(M)]
    p = [mp.mpf(0)]*(L+1)
    for n in range(L+1):
        for j in range(min(n,M)+1):
            p[n] += q[j]*c[n-j]
    return p,q

def pade_eval(p,q,z):
    num = sum(p[k]*z**k for k in range(len(p)))
    den = sum(q[k]*z**k for k in range(len(q)))
    return num/den

# Borel transform series: B(z) = sum a_k/k! z^k ; Pade [L/M]; Laplace integrate
def borel_sum(aco, tval, L, M):
    bb = [aco[k]/mp.factorial(k) for k in range(L+M+1)]
    p,q = pade(bb,L,M)
    # ∫_0^∞ e^{-u} P(t u) du  -- t may be negative; integrate in pieces
    f = lambda u: mp.e**(-u)*pade_eval(p,q,tval*u)
    return mp.quad(f, [0, 20, 60, mp.inf])

print("Padé-Borel resummation of S(1+d)/S(1) in the d-variable:")
print(f"{'d':>7} {'resummed':>25} {'direct/S1':>25} {'digits':>8}")
for dstr in ["0.001","0.003","0.01","0.03","0.1","-0.001","-0.003","-0.01","-0.03","-0.1"]:
    d = mp.mpf(dstr)
    ref = mp.mpf(direct[dstr]["N1"])/S1
    best = None
    for (L,M) in [(40,40),(50,50),(59,60),(45,45)]:
        try:
            v = borel_sum(dco, d, L, M)
        except Exception as e:
            continue
        dig = -mp.log10(abs(v-ref)/abs(ref)+mp.mpf('1e-90'))
        if best is None or dig>best[1]: best=(v,dig,(L,M))
    v,dig,LM = best
    print(f"{dstr:>7} {mp.nstr(v,22):>25} {mp.nstr(ref,22):>25} {mp.nstr(dig,5):>8}  {LM}")
