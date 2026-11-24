import sympy as sp
import mpmath as mp
import pickle
mp.mp.dps = 55
pi = mp.pi; s2 = mp.sqrt(2)

# ---- build c_n(k) for n<=10 from D~ tables ----
Dt10 = pickle.load(open("Dt10.pkl","rb"))
# reuse Dt dict structure from step_15: reconstruct programmatically
Dt = {}
exec(open("step_15.py").read().split("# verify against")[0].replace("mp.dps = 60","mp.dps=55").replace("import pickle","import pickle"))
for r, coeffs in Dt10.items():
    monos = [('z',10), ('p',(3,7)), ('p',(5,5)), ('p',(2,3,5)), ('p',(2,2,3,3))]
    Dt[(10, r)] = {mo: sp.Rational(c) for mo, c in zip(monos, coeffs)}
def monovalue(mo):
    if mo[0]=='z': return mp.zeta(mo[1])
    return mp.fprod([mp.zeta(i) for i in mo[1]])
def Q(n, uval):
    tot = mp.mpf(0)
    for r in range((n-2)//2+1):
        q = sp.Rational((-1)**r, sp.factorial(2*r))
        tot += mp.mpf(q) * (4*uval)**r * mp.fsum(mp.mpf(str(c))*monovalue(mo) for mo,c in Dt[(n,r)].items())
    return tot
def cn(n, k): return (-1)**k * s2/pi * Q(n, k*k)
def A(k): return -(-1)**k*4*s2/(pi*(16*k*k-1))

# ---- direct theta-quadrature of a_k (independent of the s-integral machinery) ----
def a_direct(k, tau, dps_check=None):
    q = 1/tau
    h = lambda t: (mp.cos(t)**q + mp.sin(t)**q)**tau  # on [0,pi/4], cos,sin >= 0
    pts = [0, pi/8, pi/4-8*tau, pi/4-2*tau, pi/4-tau/2, pi/4]
    return 8/pi*mp.quad(lambda t: h(t)*mp.cos(4*k*t), pts)

print("TASK 3: a_k by direct quadrature vs series sum_{n=2}^{10} c_n tau^n, tau=0.01 and 0.005")
for tau in [mp.mpf('0.01'), mp.mpf('0.005')]:
    print(f"--- tau = {tau} ---")
    for k in range(6):
        aq = a_direct(k, tau)
        ser = A(k) + mp.fsum(cn(n,k)*tau**n for n in range(2,11))
        diff = aq - ser
        print(f"k={k}: a_k = {mp.nstr(aq,42)}")
        print(f"     ser = {mp.nstr(ser,42)}")
        print(f"     diff = {mp.nstr(diff,4)}   (~ expected O(tau^11) = {mp.nstr(tau**11,2)} x c_11-ish)")
