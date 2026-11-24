import os
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
common = r'''"""RST1 shared numerics: exact rest and the assembled analytic bound.
M_k - M_main = -J_k + (2/pi) Im(T1 + T2 - T3)   [CX6 (8)]
  T1 = int_Gamma H K A_k ds   (dominant; q=2 pieces give e^{-sqrt(2 pi N)})
  T2 = int_Gamma H (K-E) Q_k ds
  T3 = int_0^{1/2} H E Q_k ds
Bound target: C N^{-1/2} e^{-sqrt(2 pi N)}.
"""
import mpmath as mp
from mpmath import mpf, pi, sqrt, exp, log
mp.mp.dps = 60

j_ = 1j
omega = exp(j_*pi/4)
a = pi/(2*sqrt(2))
sq2 = sqrt(2)
sqpi = sqrt(pi)
sq2pi = sqrt(2*pi)

def H(s): return mp.beta(s/2, mpf(3)/2)/2
def s_of(r): return mpf('0.5') + r*omega
def x_of(r): return mpf('0.5') + r/sq2

# ---- F_H bounds (Lemmas R2/R3): F_H(j,q) = int_0^inf e^{-a r} |H| rho^j/|s+q|^{j+1} dr
def FH_gauss(j, q):
    """Lemma R2: requires r_* >= 2 c0."""
    r_ = sqrt(2*j*q/pi); c0 = 1 + q/sq2
    assert r_ >= 2*c0, (j, q, r_, c0)
    spj = sqrt(pi*j*q)
    left = 2**mpf('1.5')/(3*a) * exp(-spj/4)
    mid  = (4*sq2)**mpf('1.5') * sqrt(125/(2*sq2)) / r_
    right= (sq2*r_)**mpf('-1.5') / a
    return sqpi/q * (left + mid + right) * exp(a*c0/2) * exp(-spj)

def FH_box(j, q):
    """Lemma R3: universal."""
    r_ = sqrt(2*j*q/pi); c0 = 1 + q/sq2
    return 2*sq2/q * (2*r_*exp(a*c0/2) + 1/a) * exp(-sqrt(pi*j*q))

# ---- T1 assembly: (2/pi)|T1| <= (8/pi)[ F_H(k,2) + sum_{i>=1} 2^i/i! F_H(k-i,i+1) ]
def T1_bound(k):
    i0 = FH_gauss(k, 2) if k >= 19 else FH_box(k, 2)
    i1 = 2*(FH_gauss(k-1, 2) if k-1 >= 19 else FH_box(k-1, 2))
    # i in [2, k/2]: box with exponent lump e^{-sqrt(3 pi (k-2))}
    S = mpf(0)
    for i in range(2, k//2+1):
        S += 2**i/mp.factorial(i)*FH_box(k-i, i+1)
    # i > k/2: box, drop exponent, CX6(16) coefficient trick done outside (here direct)
    T = mpf(0)
    for i in range(k//2+1, k+1):
        T += 2**i/mp.factorial(i)*FH_box(k-i, i+1)
    return (8/pi)*(i0 + i1 + S + T)

def T2_bound(k):
    """CX6 (17): (2/pi)|T2| <= T0 e^{-sqrt(9 pi/4) k}, T0 = 32/(3 pi a) e^{3a/2+3pi/4}."""
    T0 = 32/(3*pi*a)*exp(3*a/2 + 3*pi/4)
    return T0*exp(-sqrt(9*pi/4)*k)

def T3J_bound(k):
    """(2/pi)|T3| <= 2/(pi k) 3^{-k}  [CX6 (18)];  |J_k| <= 2^{-(N-1)}."""
    return 2/(pi*k)*3**(-k) + 2**(-k)

def rest_bound(k):
    """Analytic bound on |M_k - M_main| (raw, in k)."""
    return T1_bound(k) + T2_bound(k) + T3J_bound(k)

def target(k):
    N = k+1
    return N**mpf('-0.5')*exp(-sq2pi*sqrt(N))

# exact rest
import re
def load_mk():
    out={}
    for line in open("<project>/research/08_theory/D13f_sweep/mk_deg9.tsv"):
        p=line.rstrip("\n").split("\t")
        if not p[0].isdigit(): continue
        mid=re.match(r"\[([-0-9.eE+]+)",p[1]).group(1)
        out[int(p[0])]=mpf(mid)
    return out

def Mmain(k):
    N=k+1
    rs=mp.sqrt(2*N/pi)
    f=lambda r: exp(j_*pi*(r*omega)/2)*(r*omega/(1+r*omega))**N*H(r*omega)/(r*omega)*omega
    pts=[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]
    return 2/pi*mp.quad(f,pts).imag
'''
open('rst1_common.py','w').write(common)
step = r'''"""RST1 step_10: validate the assembled analytic bound against the exact rest."""
import mpmath as mp
from rst1_common import *
rows = load_mk()
print(" k      |rest|            bound(k)         C_eff=bound/target   ratio |rest|/bound")
for k in [99,120,150,200,250,300,350,400]:
    mm = Mmain(k); rest = abs(rows[k]-mm)
    b = rest_bound(k); t = target(k)
    print(f"{k:4d}  {mp.nstr(rest,6):>12s}   {mp.nstr(b,6):>12s}   {mp.nstr(b/t,8):>10s}   {mp.nstr(rest/b,8)}")
# constant C as max of bound/target over k>=99 (sample)
mx=0; mk=None
for k in range(99,401,2):
    r=rest_bound(k)/target(k)
    if r>mx: mx=r; mk=k
print("max bound/target over k in [99,400] step2:", mp.nstr(mx,8), "at k=",mk)
'''
open('step_10.py','w').write(step)
import subprocess
r = subprocess.run(['bash','-c','cd <project>/research/08_theory/RST1_rest_exponent && nice <project>/.venv/bin/python step_10.py'], capture_output=True, text=True, timeout=2400)
print(r.stdout)
print(r.stderr[-1500:] if r.stderr else '')
