import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import mpmath as mp
mp.mp.dps = 110
import json

PR = json.load(open("padeR9090.json"))
numR = [mp.mpf(v) for v in PR['num']]; denR = [mp.mpf(v) for v in PR['den']]
def pevalR(t):
    nv = numR[0]; dv = denR[0]
    for k in range(1, len(numR)): nv += numR[k]*t**k
    for k in range(1, len(denR)): dv += denR[k]*t**k
    return nv/dv
Hfun = lambda z: mp.sqrt(2)*mp.sqrt(1+mp.exp(-2*z))

def de_quad(phi, h, Tlo, Thi):
    tot = mp.mpf(0)
    klo = int(Tlo/h)+1; khi = int(Thi/h)+1
    for k in range(-klo, khi+1):
        t = k*h
        sh = mp.sinh(t); u = mp.exp(mp.pi/2*sh)
        up = mp.pi/2*mp.cosh(t)*u
        tot += phi(u)*up
    return tot*h

refs = {k: mp.mpf(v[0].split(' +/-')[0].strip('[]')) for k, v in json.load(open("iref_certified.json")).items()}
Tlo = mp.asinh(2/mp.pi*mp.log(mp.mpf('1e60')))
Thi = mp.asinh(2/mp.pi*mp.log(160))

print("eta=+0.3 and eta=-0.3: direct [90/90] Borel-Pade digits")
for eta in ['0.3','-0.3']:
    e = mp.mpf(eta)
    S1 = de_quad(lambda u: mp.exp(-u)*pevalR(e*u) + 0*mp.exp(-u)*Hfun(e*u), mp.mpf(1)/64, Tlo, Thi)
    # direct (full B Pade): rebuild peval from pade9090
P = json.load(open("pade9090.json"))
num = [mp.mpf(v) for v in P['num']]; den = [mp.mpf(v) for v in P['den']]
def peval(t):
    nv = num[0]; dv = den[0]
    for k in range(1, len(num)): nv += num[k]*t**k
    for k in range(1, len(den)): dv += den[k]*t**k
    return nv/dv
for eta in ['0.3','-0.3']:
    e = mp.mpf(eta)
    S = de_quad(lambda u: mp.exp(-u)*peval(e*u), mp.mpf(1)/64, Tlo, Thi)
    err = abs(S-refs[eta])
    print(f"eta={eta}: direct [90/90] digits = {mp.nstr(-mp.log10(err),6)}")
'''
with open(wd+"/step_45.py","w") as fp: fp.write(code)
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_45.py"], capture_output=True, text=True, cwd=wd, timeout=900)
print(r.stdout[-1500:]); print("STDERR:", r.stderr[-500:]); print("elapsed", time.time()-t0)
