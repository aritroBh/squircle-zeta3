
import mpmath as mp
mp.mp.dps = 160
import json

def de_quad(phi, h, Tlo, Thi):
    tot = mp.mpf(0)
    klo = int(Tlo/h)+1; khi = int(Thi/h)+1
    for k in range(-klo, khi+1):
        t = k*h
        sh = mp.sinh(t); u = mp.exp(mp.pi/2*sh)
        up = mp.pi/2*mp.cosh(t)*u
        tot += phi(u)*up
    return tot*h

def phi_u(eta):
    def phi(u):
        return mp.exp(-u)*(2*mp.exp(-eta*u/2)*mp.sqrt(mp.cosh(eta*u)))*(1+mp.exp(-u))**(-2-eta)
    return phi

Tlo = mp.asinh(2/mp.pi*mp.log(mp.mpf('1e80')))
Thi = mp.asinh(2/mp.pi*mp.log(200))
refs = {}
for eta in ['0.05','0.1','0.2','0.3','-0.05','-0.1','-0.2','-0.3']:
    e = mp.mpf(eta); phi = phi_u(e)
    v1 = de_quad(phi, mp.mpf(1)/64, Tlo, Thi)
    v2 = de_quad(phi, mp.mpf(1)/96, Tlo, Thi)
    refs[eta] = v2
    print(f"eta={eta}: 1/64 vs 1/96 agree: {mp.nstr(abs(v1-v2),3)}")
    print("   ", mp.nstr(v2, 80))
json.dump({k: mp.nstr(v, 130) for k, v in refs.items()}, open("iref_de160.json","w"))
