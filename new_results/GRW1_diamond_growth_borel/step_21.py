
import mpmath as mp
mp.mp.dps = 120
import json

PR = json.load(open("padeR9090.json"))
numR = [mp.mpf(v) for v in PR['num']]; denR = [mp.mpf(v) for v in PR['den']]
def pevalR(t):
    nv = numR[0]; dv = denR[0]
    for k in range(1, len(numR)): nv += numR[k]*t**k
    for k in range(1, len(denR)): dv += denR[k]*t**k
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
Tlo = mp.asinh(2/mp.pi*mp.log(mp.mpf('1e80')))
Thi = mp.asinh(2/mp.pi*mp.log(200))

Hfun = lambda z: mp.sqrt(2)*mp.sqrt(1+mp.exp(-2*z))
print("=== S(eta) = S_H + S_R  (S_H = int e^{-u} H(eta u) du) ===")
table = {}
for eta in ['0.05','0.1','0.2','0.3','-0.05','-0.1','-0.2','-0.3']:
    e = mp.mpf(eta)
    SH1 = de_quad(lambda u: mp.exp(-u)*Hfun(e*u), mp.mpf(1)/48, Tlo, Thi)
    SH2 = de_quad(lambda u: mp.exp(-u)*Hfun(e*u), mp.mpf(1)/64, Tlo, Thi)
    SR1 = de_quad(lambda u: mp.exp(-u)*pevalR(e*u), mp.mpf(1)/48, Tlo, Thi)
    SR2 = de_quad(lambda u: mp.exp(-u)*pevalR(e*u), mp.mpf(1)/64, Tlo, Thi)
    S = SH2 + SR2
    err = abs(S - refs[eta])
    digits = -mp.log10(err)
    table[eta] = (mp.nstr(S, 80), mp.nstr(digits, 5))
    print(f"eta={eta}: SH-agree={mp.nstr(abs(SH1-SH2),3)} SR-agree={mp.nstr(abs(SR1-SR2),3)}")
    print(f"   S={mp.nstr(S,45)}")
    print(f"   |S-I|={mp.nstr(err,3)}  digits={mp.nstr(digits,5)}  (paper target: 44/31/25/30 at 0.1/0.2/0.3/-0.2)")
json.dump(table, open("pb_split_sums.json","w"))
