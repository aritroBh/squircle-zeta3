import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import mpmath as mp
mp.mp.dps = 350
import json

d = json.load(open("<project>/research/08_theory/D11b2_stokes_jump/step2_num450_N180.json"))
f = [mp.mpf(v) for v in d['f450']]
b = [f[n]/mp.factorial(n) for n in range(len(f))]
def pade_coeffs(a, L, M):
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k = L + i - j
            A[i, j] = a[k] if k >= 0 else mp.mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = (A**-1)*rhs
    den = [mp.mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j]*a[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
    return num, den

refs = {k: mp.mpf(v[0].split(' +/-')[0].strip('[]')) for k, v in json.load(open("iref_certified.json")).items()}

def de_quad(phi, h, Tlo, Thi):
    tot = mp.mpf(0)
    klo = int(Tlo/h)+1; khi = int(Thi/h)+1
    for k in range(-klo, khi+1):
        t = k*h
        sh = mp.sinh(t); u = mp.exp(mp.pi/2*sh)
        up = mp.pi/2*mp.cosh(t)*u
        tot += phi(u)*up
    return tot*h
Tlo = mp.asinh(2/mp.pi*mp.log(mp.mpf('1e60')))
Thi = mp.asinh(2/mp.pi*mp.log(160))

import sys
for M in [90, 85, 80, 75, 70, 60]:
    num, den = pade_coeffs(b, M, M)
    numS = [mp.nstr(v, 340) for v in num]; denS = [mp.nstr(v, 340) for v in den]
    mp.mp.dps = 110
    num = [mp.mpf(v) for v in numS]; den = [mp.mpf(v) for v in denS]
    def peval(t):
        nv = num[0]; dv = den[0]
        for k in range(1, len(num)): nv += num[k]*t**k
        for k in range(1, len(den)): dv += den[k]*t**k
        return nv/dv
    for eta in ['-0.2', '0.1']:
        e = mp.mpf(eta)
        # DE quadrature
        S1 = de_quad(lambda u: mp.exp(-u)*peval(e*u), mp.mpf(1)/64, Tlo, Thi)
        # adaptive quad
        S2 = mp.quad(lambda t: mp.exp(-t)*peval(e*t), [0, 2, 8, 20, 60])
        ref = refs[eta]
        print(f"[{M}/{M}] eta={eta}: DE: {mp.nstr(-mp.log10(abs(S1-ref)),6)} digits, adaptive: {mp.nstr(-mp.log10(abs(S2-ref)),6)} digits  |S1-S2|={mp.nstr(abs(S1-S2),2)}")
    mp.mp.dps = 350
'''
with open(wd+"/step_44.py","w") as fp: fp.write(code)
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_44.py"], capture_output=True, text=True, cwd=wd, timeout=1700)
print(r.stdout[-2500:]); print("STDERR:", r.stderr[-800:]); print("elapsed", time.time()-t0)
