
import mpmath as mp
mp.mp.dps = 110
import json

P = json.load(open("pade9090.json"))
num = [mp.mpf(v) for v in P['num']]; den = [mp.mpf(v) for v in P['den']]
def peval(t):
    nv = num[0]; dv = den[0]
    for k in range(1, len(num)): nv += num[k]*t**k
    for k in range(1, len(den)): dv += den[k]*t**k
    return nv/dv

def de_quad(phi, h, Tlo, Thi):
    tot = mp.mpf(0)
    klo = int(Tlo/h)+1; khi = int(Thi/h)+1
    for k in range(-klo, khi+1):
        t = k*h
        sh = mp.sinh(t); u = mp.exp(mp.pi/2*sh)
        up = mp.pi/2*mp.cosh(t)*u
        tot += phi(u)*up
    return tot*h

refs = {k: mp.mpf(v) for k, v in json.load(open("iref_de160.json")).items()}

# pole scan: PB denominator on the real axis (min |den|) and near positive axis
import numpy as np
den_d = np.array([float(v) for v in den[::-1]])  # highest degree first
rts = np.roots(den_d)
rts_sorted = sorted(rts, key=abs)
print("smallest |poles| (double precision, indicative):")
for z in rts_sorted[:8]:
    print(f"   {z:.6f}   |z|={abs(z):.6f}   (pi/2={float(mp.pi/2):.6f})")
minreal = min(abs(peval_den) if False else abs(sum(den[k]*x**k for k in range(len(den)))) for x in [i*0.5 for i in range(-80, 81)])
print("min |den(x)| on real x in [-40,40]:", mp.nstr(minreal, 6))

print()
print("=== Pade-Borel sum S(eta) = int_0^inf e^{-u} PB(eta u) du vs reference ===")
Tlo = mp.asinh(2/mp.pi*mp.log(mp.mpf('1e60')))
Thi = mp.asinh(2/mp.pi*mp.log(160))
res = {}
for eta in ['0.05','0.1','0.2','0.3','-0.05','-0.1','-0.2','-0.3']:
    e = mp.mpf(eta)
    phi = lambda u: mp.exp(-u)*peval(e*u)
    v1 = de_quad(phi, mp.mpf(1)/48, Tlo, Thi)
    v2 = de_quad(phi, mp.mpf(1)/64, Tlo, Thi)
    err = abs(v2 - refs[eta])
    res[eta] = mp.nstr(v2, 90)
    digits = -mp.log10(err) if err > 0 else mp.mpf(99)
    print(f"eta={eta}: S={mp.nstr(v2, 40)}")
    print(f"         quad-consistency={mp.nstr(abs(v1-v2),3)}  |S-I|={mp.nstr(err,3)}  digits={mp.nstr(digits,5)}")
json.dump(res, open("pb_sums.json","w"))
