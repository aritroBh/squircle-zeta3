lllk = r'''"""lllkernel.py -- approximate integer kernel of a real matrix via LLL (python-flint).
Given A (M x N) mpf rows, K scale: returns short integer vectors q with A q ~ 0."""
import mpmath as mp
from flint import fmpz_mat, fmpz

def int_kernel(A, Kexp=45):
    """A: list of rows (each a list of mpf), M x N. Returns list of integer vectors q (len N)."""
    M = len(A); N = len(A[0])
    K = mp.power(10, Kexp)
    # B = [I_N | round(K * A^T)] : N rows, N+M cols
    B = fmpz_mat(N, N + M)
    for i in range(N):
        B[i, i] = fmpz(1)
    for j in range(M):
        for i in range(N):
            x = mp.nint(K*A[j][i])
            B[i, N + j] = fmpz(int(x))
    BL = B.lll()
    out = []
    for i in range(N):
        q = [int(BL[i, j]) for j in range(N)]
        tail = [int(BL[i, N + j]) for j in range(M)]
        out.append((q, max(abs(t) for t in tail) if tail else 0))
    return out
'''
open('lllkernel.py','w').write(lllk)
step25 = r'''"""step_25.py -- LLL integer-kernel identification: test on weight two (recover T), then weight-three slope."""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 80
from lllkernel import int_kernel

def T_alpha_c(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def Rw2(w): return 1/(2*w) + 1/(2*(w-1)) - 1/(w+1) - 1/(2*(w+r))
    return mp.quad(lambda w: Rw2(w)*lt(w), [r, mp.sqrt(r), 1])

def Jval(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def lmt(w): return mp.log((1-w)*(w+r)/(2*(1-r)*w))
    def f(w):
        return (1/(w+1) - 1/(w+r))*lt(w)**2 + (1/w - 2/(w+1))*lt(w)*lmt(w)
    return mp.quad(f, [r, mp.sqrt(r), 1])

# ---- test: weight two ----
angles = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7, mp.mpf(5)/6]
def basis2(c):
    return [('pi^2', mp.pi**2), ('Li2(c)', mp.polylog(2, c)), ('Li2(-c)', mp.polylog(2, -c)), ('log2^2', mp.log(2)**2)]
A = [[T_alpha_c(c)] + [x for _, x in basis2(c)] for c in angles]
names = ['T'] + [n for n, _ in basis2(angles[0])]
for q, tail in int_kernel(A, Kexp=45)[:3]:
    if tail < 10**10 and any(q):
        print("W2 relation:", [(cf, names[i]) for i, cf in enumerate(q) if cf], " tail:", tail)

# ---- weight three slope ----
def basis3(c):
    L = mp.log; A3 = lambda x: mp.polylog(3, x)
    return [('zeta3', mp.zeta(3)),
            ('Li3(c)', A3(c)), ('Li3(-c)', A3(-c)), ('Li3(1-c)', A3(1-c)), ('Li3(1/(1+c))', A3(1/(1+c))),
            ('Li3((1+c)/2)', A3((1+c)/2)), ('Li3((1-c)/(1+c))', A3((1-c)/(1+c))), ('Li3(-(1-c)/(1+c))', A3(-(1-c)/(1+c))),
            ('Li3(2c/(1+c))', A3(2*c/(1+c))),
            ('pi^2 log2', mp.pi**2*L(2)), ('pi^2 logc', mp.pi**2*L(c)), ('pi^2 log(1+c)', mp.pi**2*L(1+c)), ('pi^2 log(1-c)', mp.pi**2*L(1-c)),
            ('log2^3', L(2)**3), ('logc^3', L(c)**3), ('log(1+c)^3', L(1+c)**3), ('log(1-c)^3', L(1-c)**3),
            ('log2^2 logc', L(2)**2*L(c)), ('log2 logc^2', L(2)*L(c)**2), ('logc^2 log(1+c)', L(c)**2*L(1+c)),
            ('logc log(1+c)^2', L(c)*L(1+c)**2), ('log(1+c)^2 log(1-c)', L(1+c)**2*L(1-c)), ('log2 logc log(1+c)', L(2)*L(c)*L(1+c))]

print("evaluating slope...")
def slope(c): return (1-c**2)/(2*c**2)*Jval(c)
angles3 = [mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3, mp.mpf(4)/5, mp.mpf(2)/7, mp.mpf(5)/6, mp.mpf(1)/5, mp.mpf(3)/7]
A3m = [[slope(c)] + [x for _, x in basis3(c)] for c in angles3]
names3 = ['slope'] + [n for n, _ in basis3(angles3[0])]
print("relations:")
for q, tail in int_kernel(A3m, Kexp=45)[:4]:
    if tail < 10**12 and any(q):
        print("W3 relation:", [(cf, names3[i]) for i, cf in enumerate(q) if cf], " tail:", tail)
print("STEP25 DONE")
'''
open('step_25.py','w').write(step25)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_25.py'], capture_output=True, text=True, timeout=1700)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else "RC=0")
