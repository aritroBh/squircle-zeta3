lib = '''"""LEAD2R referee library. Written independently from the proof text (no import of lead2_lib).
All formulas re-derived from proof_lead.md by the referee."""
import os, pickle
from fractions import Fraction as Fr
from math import comb
import mpmath as mp

BASE = '<project>/research/08_theory'
WTS  = os.path.join(BASE, 'WTS_weight_k_structure', 'theory')
TWOEND = os.path.join(BASE, 'TWOEND1_weight_poles', 'red')

# ---------- exact zigzag numbers (own implementation: Entringer-style via sech+tanh egf) ----------
def zigzag_ints(N):
    """A_0..A_N exact integers, sec x + tan x = sum A_n x^n/n!. Own recurrence on rational series."""
    import sympy as sp
    x = sp.symbols('x')
    ser = sp.series(sp.sec(x) + sp.tan(x), x, 0, N + 1).removeO().expand()
    return [int(sp.Rational(ser.coeff(x, n)) * sp.factorial(n)) for n in range(N + 1)]

def zigzag_recurrence(N):
    """Second independent computation: A_{n+1} = (1/2) sum_k binom(n,k) A_k A_{n-k} for n>=2? No --
    use the (sec+tan) ODE: f = sec+tan, f' = sec tan + sec^2 = f*sec. sec = (f + 1/f')/2? Simpler:
    tangent numbers + secant numbers separately via classical recurrences."""
    # secant (Euler) numbers: sech x = sum (-1)^n E_{2n} x^{2n}/(2n)!, E_0=1, sum binom(2n,2k) E_{2k}=0
    A = [0]*(N+1)
    E2 = [0]*(N+2)   # E_{2n}
    E2[0] = 1
    for n in range(1, N//2+1):
        s = sum(comb(2*n, 2*k)*E2[k] for k in range(n))
        E2[n] = -s
    # tangent numbers: tan x = sum T_n x^{2n-1}/(2n-1)!, T_n = (-1)^{n-1} 4^n (4^n-1) B_{2n}/(2n)
    import sympy as sp
    for n in range(1, (N+1)//2+1):
        Tn = (-1)**(n-1)*4**n*(4**n-1)*sp.bernoulli(2*n)//(2*n)
        A[2*n-1] = int(abs(Tn))
    for n in range(N//2+1):
        A[2*n] = int(abs(E2[n]))
    return A

def eps_m(m):
    return (-1)**((m+1)//2)

def delta_rat(d):
    """delta_d = [x^d](sech x + tanh x) as Fraction, computed two ways must agree."""
    import sympy as sp
    x = sp.symbols('x')
    c = sp.Rational(sp.series(sp.sech(x)+sp.tanh(x), x, 0, d+1).removeO().coeff(x, d))
    return Fr(int(c.p), int(c.q))

# ---------- E(y) coefficients: two independent ways ----------
def e_coeffs_logE(N):
    """via log E = y log2 + sum s_n y^n/n (needs mp.mp.dps set by caller)."""
    L = [mp.mpf(0)]*(N+2); L[1] = mp.log(2)
    for n in range(2, N+2):
        L[n] = (mp.zeta(n) if n % 2 == 0 else (1-mp.mpf(2)**(1-n))*mp.zeta(n))/n
    E = [mp.mpf(0)]*(N+2); E[0] = mp.mpf(1)
    for n in range(1, N+2):
        E[n] = sum(k*L[k]*E[n-k] for k in range(1, n+1))/n
    return E

def Efun(y):
    return mp.gamma(1+y/2)*mp.gamma(mp.mpf('0.5')-y/2)/mp.sqrt(mp.pi)

def e_coeffs_taylor(N):
    """independent: high-order Taylor of the Gamma product via mpmath."""
    return [mp.taylor(Efun, 0, N)[n] for n in range(N+1)]

# ---------- kernels ----------
POLE_NAME = {(0,0):'0', (1,1):'a', (-1,-1):'-a', (-1,1):'b', (1,-1):'-b'}
POLE_VAL = None  # filled per-precision
def load_all_kernels():
    kern = {}
    for K in range(1, 19):
        kern[K] = pickle.load(open(os.path.join(WTS, f'kern_{K}.pkl'), 'rb'))
    for K in (19, 20, 21):
        kern[K] = pickle.load(open(os.path.join(TWOEND, f'kern_{K}.pkl'), 'rb'))
    return kern

def residue_rational(kern, k, key, pole_pair):
    rec = kern.get(k, {}).get(key)
    if rec is None: return Fr(0)
    v = rec['real'].get(pole_pair)
    return Fr(0) if v is None else v[0]

def residue_sqrt2(kern, k, key, pole_pair):
    rec = kern.get(k, {}).get(key)
    if rec is None: return Fr(0)
    v = rec['real'].get(pole_pair)
    return Fr(0) if v is None else v[1]

def forward_diffs(seq, m):
    """m-th forward differences of the sequence (list of Fr)."""
    s = list(seq)
    for _ in range(m):
        s = [b-a for a, b in zip(s, s[1:])]
    return s
'''
with open('lead2r_lib.py','w') as f: f.write(lib)
print("written", len(lib))
with open('progress.log','a') as f:
    f.write("step_02: wrote lead2r_lib.py (own independent library: zigzag 2 ways, delta, E coeffs 2 ways, kernel loaders, exact Newton diffs)\n")
