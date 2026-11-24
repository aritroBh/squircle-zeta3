
import mpmath as mp
mp.mp.dps = 120
import json

# ---- manual double-exponential (tanh-sinh) quadrature on [0,inf): u = exp(pi/2 sinh t) ----
def de_quad(phi, h, T):
    tot = mp.mpf(0)
    kmax = int(T/h)+1
    for k in range(-kmax, kmax+1):
        t = k*h
        sh = mp.sinh(t); u = mp.exp(mp.pi/2*sh)
        up = mp.pi/2*mp.cosh(t)*u
        tot += phi(u)*up
    return tot*h

def phi_u(eta):
    def phi(u):
        return mp.exp(-u)*(2*mp.exp(-eta*u/2)*mp.sqrt(mp.cosh(eta*u)))*(1+mp.exp(-u))**(-2-eta)
    return phi

print("=== references via manual tanh-sinh, two meshes ===")
refs = {}
for eta in ['0.05','0.1','0.2','0.3','-0.05','-0.1','-0.2','-0.3']:
    e = mp.mpf(eta)
    phi = phi_u(e)
    # tail beyond u=U: e^{-U}; choose U=140: e^{-140} ~ 1e-61
    T = mp.asinh(2/mp.pi*mp.log(140))
    v1 = de_quad(phi, mp.mpf(1)/32, T)
    v2 = de_quad(phi, mp.mpf(1)/64, T)
    refs[eta] = v2
    print(f"eta={eta}: h=1/32 vs 1/64 agree {mp.nstr(abs(v1-v2),3)}")
    print("   ", mp.nstr(v2, 70))
json.dump({k: mp.nstr(v, 100) for k, v in refs.items()}, open("iref_de.json","w"))
