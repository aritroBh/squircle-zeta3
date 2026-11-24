import sys, json, re, itertools
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor4_sym import *
import sympy as sp

h = Fr(1,2)
def cmul2(a, b): return (a[0]*b[0] + 2*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def cadd(a, b): return (a[0]+b[0], a[1]+b[1])
def cneg(a): return (-a[0], -a[1])
ONE=(Fr(1),Fr(0)); ZERO=(Fr(0),Fr(0)); SQ2=(Fr(0),Fr(1)); A1=(Fr(1),Fr(1))  # 1, 0, sqrt2, 1+sqrt2
def cdiv(a, b):
    # a/b in Q(sqrt2): b=(p,q), 1/b = (p,-q)/(p^2-2q^2)
    p, q = b; n = p*p - 2*q*q
    return cmul2(a, (p/n, -q/n))
def logvec2(z):
    """z in Q(sqrt2), return (j,m) with z = 2^{j/2} (1+sqrt2)^m as exponent vector over (2,1+sqrt2):
    i.e. z = 2^e0 * (1+sqrt2)^e1 with e0=j/2. Verify exactly; else raise."""
    if z == ZERO: return None
    neg = False
    zz = z
    # sign: use float for sign only
    import math
    if float(z[0]) + float(z[1])*math.sqrt(2) < 0:
        neg = True; zz = cneg(zz)
    # divide out 1+sqrt2 powers: (1+sqrt2)^m
    m = 0
    while True:
        f = float(zz[0]) + float(zz[1])*math.sqrt(2)
        if f > 1.5: zz = cdiv(zz, A1); m += 1
        elif f < 0.6: zz = cmul2(zz, A1); m -= 1
        else: break
    # now zz should be a power of 2 times maybe residual; divide by sqrt2^k via norm
    N = cmul2(zz, (zz[0], -zz[1]))[0]
    if N < 0: N = -N
    # N must be a power of 2
    import math
    if N == 0: raise ValueError('zero norm')
    jf = math.log2(float(N))
    if abs(jf - round(jf)) > 1e-9: raise ValueError(f'not S-unit: {z}, norm {N}')
    j2 = int(round(jf))  # norm exponent: N = 2^{j2}
    # zz/2^{j2/2}: zz = 2^a * unit with unit norm +-1 and in <0.6,1.5>: unit = +-1
    # a = j2/2
    if j2 % 2 == 0:
        cand = (Fr(2)**(j2//2), Fr(0))
        ok = (zz == cand)
    else:
        # zz = sqrt2 * 2^((j2-1)/2) * +-1 ?
        cand = cmul2(SQ2, (Fr(2)**((j2-1)//2), Fr(0)))
        ok = (zz == cand)
    if not ok:
        raise ValueError(f'not pure S-unit: {z} -> zz={zz}, m={m}, N={N}')
    # exponent vector over (2, 1+sqrt2): 2 = 2^1; sqrt2 = 2^{1/2}
    e0 = Fr(j2, 2)
    return (e0, Fr(-m))  # zz = 2^{e0} (1+sqrt2)^{-m}? sign check: we divided by A1 m times -> zz*A1^{...}: z = zz * A1^m... 

# verify logvec2 on known elements
for z, expect in [((h, Fr(0)), None), ((h, Fr(-1)), None), ((Fr(0), Fr(-1)), None)]:
    pass

# Simpler robust approach: solve z = 2^{e0}(1+sqrt2)^{e1} via norm + one embedding
def sunit_vec(z):
    """z in Q(sqrt2)^*. Return (e0,e1) with |z| = 2^e0 (1+sqrt2)^e1 exactly, or None."""
    import math
    if z == ZERO: return None
    sgn = float(z[0]) + float(z[1])*math.sqrt(2)
    zz = z if sgn > 0 else cneg(z)
    N = cmul2(zz, (zz[0], -zz[1]))[0]
    if N <= 0: return None
    jf = math.log2(float(N))
    if abs(jf - round(jf)) > 1e-9: return None
    j2 = int(round(jf))
    # |z| = 2^{j2/2} u with u unit of norm +1: u = (1+sqrt2)^m
    # m from archimedean: |z| = 2^{j2/2} (1+sqrt2)^m
    lf = (math.log(abs(float(zz[0])+float(zz[1])*math.sqrt(2))) - (j2/2)*math.log(2)) / math.log(1+math.sqrt(2))
    if abs(lf - round(lf)) > 1e-9: return None
    m = int(round(lf))
    # verify: zz ?= 2^{j2/2} (1+sqrt2)^m  (with 2^{1/2} = sqrt2 when j2 odd)
    base2 = cmul2(SQ2, (Fr(2)**((j2-1)//2), Fr(0))) if j2 % 2 else (Fr(2)**(j2//2), Fr(0))
    if j2 < 0:
        base2 = cdiv((Fr(1),Fr(0)), cmul2(SQ2, (Fr(2)**((-j2-1)//2), Fr(0)))) if (-j2) % 2 else (Fr(1)/Fr(2)**((-j2)//2), Fr(0))
    u = (Fr(1), Fr(0))
    for _ in range(abs(m)):
        u = cmul2(u, A1) if m > 0 else cdiv(u, A1)
    cand = cmul2(base2, u)
    if cand != zz: return None
    return (Fr(j2,2), Fr(m))

# tests
print(sunit_vec((h, Fr(0))), "expect (-1/2, 0) for 1/sqrt2")
print(sunit_vec((h, Fr(-1))), "expect (-1/2,-1) for 1-1/sqrt2")
print(sunit_vec((h, Fr(1))), "expect (-1/2, 1) for 1+1/sqrt2")
print(sunit_vec((Fr(3), Fr(0))), "expect None (3)")
