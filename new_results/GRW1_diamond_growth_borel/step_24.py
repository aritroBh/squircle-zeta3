
import mpmath as mp
mp.mp.dps = 70
import json

PR = json.load(open("padeR9090.json"))
numR = [mp.mpf(v) for v in PR['num']]; denR = [mp.mpf(v) for v in PR['den']]
def pevalR(t):
    nv = numR[0]; dv = denR[0]
    for k in range(1, len(numR)): nv += numR[k]*t**k
    for k in range(1, len(denR)): dv += denR[k]*t**k
    return nv/dv
Hfun = lambda z: mp.sqrt(2)*mp.sqrt(1+mp.exp(-2*z))
Bfun = lambda z: Hfun(z) + pevalR(z)

alpha = mp.pi/6
def lateral(eps, sign):
    # S_{sign*alpha}(i eps) = (e^{sign i alpha}/eps) int_0^R e^{-(r/eps) e^{sign i alpha}} B(r e^{i(pi/2+sign alpha)}) dr
    iu = mp.mpc(0,1)
    ray = lambda r: r*mp.e**(iu*(mp.pi/2 + sign*alpha))
    w = lambda r: mp.e**(-(r/eps)*mp.e**(sign*iu*alpha))
    R = 30*eps
    pts = [0, eps, 3*eps, 8*eps, 16*eps, R]
    v = mp.quad(lambda r: w(r)*Bfun(ray(r)), pts)
    return mp.e**(sign*iu*alpha)/eps * v

def I_complex(eta):
    # I(eta) = int_0^inf e^{-u} H(eta u) K(eta,u) du  (valid any complex eta, |arg safe|)
    f = lambda u: mp.exp(-u)*Hfun(eta*u)*(1+mp.exp(-u))**(-2-eta)
    return mp.quad(f, [0, 1, 4, 12, 40, 120])

print("eps, |J|, |Jpred|, ratio |.|, phase(J), phase(pred), and S_-alpha vs I(i eps)")
print("pred J = 2 sqrt(pi) e^{-3i pi/4} eps^{1/2} e^{-pi/(2 eps)}")
for eps in [mp.mpf('0.2'), mp.mpf('0.15'), mp.mpf('0.1'), mp.mpf('0.075'), mp.mpf('0.05')]:
    Sp = lateral(eps, +1)
    Sm = lateral(eps, -1)
    J = Sp - Sm
    Jpred = 2*mp.sqrt(mp.pi)*mp.e**(-3*mp.mpc(0,1)*mp.pi/4)*mp.sqrt(eps)*mp.e**(-mp.pi/(2*eps))
    Iv = I_complex(mp.mpc(0,1)*eps)
    print(f"eps={mp.nstr(eps,6)}:")
    print(f"   |J|={mp.nstr(abs(J),10)}  |Jpred|={mp.nstr(abs(Jpred),10)}  ratio={mp.nstr(abs(J)/abs(Jpred),10)}")
    print(f"   arg(J)={mp.nstr(mp.arg(J),10)}  arg(pred)={mp.nstr(mp.arg(Jpred),10)}  (-3pi/4={mp.nstr(-3*mp.pi/4,8)})")
    print(f"   S_-(i eps)={mp.nstr(Sm,18)}")
    print(f"   I(i eps) ={mp.nstr(Iv,18)}   |S_- - I|={mp.nstr(abs(Sm-Iv),3)}")
