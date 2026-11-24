
import mpmath as mp
mp.mp.dps = 100
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

def I_direct(eta):
    f = lambda u: mp.exp(-u)*Hfun(eta*u)*(1+mp.exp(-u))**(-2-eta)
    return mp.quad(f, [0, 1, 4, 12, 40, 150])

Tlo = mp.asinh(2/mp.pi*mp.log(mp.mpf('1e60')))
Thi = mp.asinh(2/mp.pi*mp.log(250))
print("eta, S_H+S_R (Borel sum), I (direct), digits")
for eta in ['0.5','0.7','0.9','-0.5','-0.7','-0.9']:
    e = mp.mpf(eta)
    SH = de_quad(lambda u: mp.exp(-u)*Hfun(e*u), mp.mpf(1)/64, Tlo, Thi)
    SR = de_quad(lambda u: mp.exp(-u)*pevalR(e*u), mp.mpf(1)/64, Tlo, Thi)
    S = SH + SR
    Iv = I_direct(e)
    err = abs(S-Iv)
    print(f"eta={eta}: S={mp.nstr(S,30)}")
    print(f"         I={mp.nstr(Iv,30)}   |S-I|={mp.nstr(err,3)}  digits={mp.nstr(-mp.log10(err),5)}")
