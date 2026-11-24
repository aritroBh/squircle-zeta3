
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
step03 = r'''
"""
SRF3 step_03: numerical Borel plane. Padé approximants of the Borel transform
B(z) = sum a_k z^k/k! (tau-series of S/S(1), K=80, 150-digit coeffs).
Poles should line up along branch cuts starting at the nearest singularities.
We check: z0=(-log2+-i pi)/2 only, or also (log2+-i pi)/2 (vertex candidates)?
"""
import mpmath as mp, json, os
mp.mp.dps = 80
base = "<project>"
h1 = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json")))
tau = [mp.mpf(x) for x in h1['eta']]
K = 80
z0e = (-mp.log(2)+mp.pi*1j)/2
z0v = ( mp.log(2)+mp.pi*1j)/2

bb = [tau[k]/mp.factorial(k) for k in range(K+1)]
# Padé [L/M] via solving the Toeplitz system (mpmath has pade? implement)
def pade_poles(c, L, M):
    # c[0..L+M]; solve for denominator q_0..q_M, q_0=1
    A = mp.matrix(M, M)
    r = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            A[i,j] = c[L-M+i+j+1]
        r[i] = -c[L+i+1]
    q = mp.lu_solve(A, r)
    q = [mp.mpf(1)] + [q[i] for i in range(M)]
    # roots of q(w)=0 with w = z: polynomial in z
    # q(z) = sum q_j z^j
    coeffs = q
    roots = mp.polyroots(coeffs[::-1])  # highest first? polyroots wants ascending
    return roots

# use mpmath's pade if available
try:
    from mpmath import pade as mp_pade
    have=True
except Exception as e:
    have=False
print("mpmath pade available:", have)

def pade_coeffs(c, L, M):
    # denominator q_0..q_M with q_M=1? standard: solve with q_0=1
    A = mp.matrix(M, M); r = mp.matrix(M,1)
    for i in range(M):
        for j in range(M):
            A[i,j] = c[L-M+i+j+1]
        r[i] = -c[L+i+1]
    qsol = mp.lu_solve(A, r)
    q = [mp.mpf(1)] + [qsol[i] for i in range(M)]
    p = [mp.mpf(0)]*(L+1)
    for n in range(L+1):
        s = mp.mpf(0)
        for j in range(min(n,M)+1):
            s += q[j]*c[n-j]
        p[n]=s
    return p,q

for L,M in [(39,40),(35,40),(30,30)]:
    p,q = pade_coeffs(bb, L, M)
    roots = mp.polyroots(list(reversed(q)))
    # polyroots expects ascending order: q[0]+q[1]z+...
    roots = mp.polyroots(q)
    dist_e = min(abs(r-z0e) for r in roots)
    dist_eb = min(abs(r-mp.conj(z0e)) for r in roots)
    dist_v = min(abs(r-z0v) for r in roots)
    dist_vb = min(abs(r-mp.conj(z0v)) for r in roots)
    print(f"Pade [{L}/{M}]: dist to edge z0+={mp.nstr(dist_e,3)}, z0-={mp.nstr(dist_eb,3)}; "
          f"to vertex z0+={mp.nstr(dist_v,3)}, z0-={mp.nstr(dist_vb,3)}")
    # print the 10 poles of smallest modulus
    roots_sorted = sorted(roots, key=abs)
    for r_ in roots_sorted[:10]:
        print("    pole at", mp.nstr(r_,12), " |z|=", mp.nstr(abs(r_),8))
'''
open(wd+"/step_03.py","w").write(step03)
import subprocess
PY = "<project>/.venv/bin/python"
r = subprocess.run(["nice","-n","10",PY,"step_03.py"],capture_output=True,text=True,cwd=wd,timeout=900)
print(r.stdout); print(r.stderr[-1500:])
