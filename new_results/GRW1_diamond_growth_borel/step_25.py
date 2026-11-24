
import mpmath as mp
mp.mp.dps = 60
import json

PR = json.load(open("padeR9090.json"))
numR = [mp.mpf(v) for v in PR['num']]; denR = [mp.mpf(v) for v in PR['den']]
def pevalR(t):
    nv = numR[0]; dv = denR[0]
    for k in range(1, len(numR)): nv += numR[k]*t**k
    for k in range(1, len(denR)): dv += denR[k]*t**k
    return nv/dv
Hfun = lambda z: mp.sqrt(2)*mp.sqrt(1+mp.exp(-2*z))
iu = mp.mpc(0,1)
alpha = mp.pi/6

def de_ray(sign, eps, fun):
    # int_0^R e^{-(r/eps) e^{sign i alpha}} fun(r e^{i(pi/2+sign alpha)}) dr, DE map on r
    Thi = mp.asinh(2/mp.pi*mp.log(40*eps))  # r up to 40 eps: e^{-40 cos alpha} ~ 9e-16
    h = mp.mpf(1)/48
    tot = mp.mpc(0)
    for k in range(-int(3.0/h)-1, int(Thi/h)+2):
        t = k*h
        sh = mp.sinh(t); r = mp.exp(mp.pi/2*sh)
        rp = mp.pi/2*mp.cosh(t)*r
        z = r*mp.e**(iu*(mp.pi/2+sign*alpha))
        w = mp.e**(-(r/eps)*mp.e**(sign*iu*alpha))
        tot += w*fun(z)*rp
    return tot*h

def lateral_R(sign, eps):
    return mp.e**(sign*iu*alpha)/eps * de_ray(sign, eps, pevalR)
def lateral_H(sign, eps):
    return mp.e**(sign*iu*alpha)/eps * de_ray(sign, eps, Hfun)

def loop_H(eps):
    # -4 i int_0^pi e^{-(pi/2+s)/eps} sqrt(sin s) e^{i(pi/4 - s/2)} ds   [k=0 cut contribution]
    f = lambda s: mp.e**(-(mp.pi/2+s)/eps)*mp.sqrt(mp.sin(s))*mp.e**(iu*(mp.pi/4 - s/2))
    v = mp.quad(f, [0, eps, 4*eps, mp.pi/2, mp.pi])
    return -4*iu*v

def I_complex(eta):
    f = lambda u: mp.exp(-u)*Hfun(eta*u)*(1+mp.exp(-u))**(-2-eta)
    return mp.quad(f, [0, 1, 4, 12, 40, 120])

print("Stokes jump across arg eta = pi/2 at eta = i eps:")
print("pred: J = 2 sqrt(pi) e^{-3 i pi/4} eps^{1/2} e^{-pi/(2 eps)}")
res = {}
for eps in [mp.mpf('0.2'), mp.mpf('0.15'), mp.mpf('0.1'), mp.mpf('0.075')]:
    SRp = lateral_R(+1, eps); SRm = lateral_R(-1, eps)
    SHm = lateral_H(-1, eps)
    lp = loop_H(eps)
    Sm = SHm + SRm            # right lateral sum (should equal I(i eps))
    Sp = Sm + lp + (SRp - SRm) # left lateral sum = I + J
    J = Sp - Sm
    Jpred = 2*mp.sqrt(mp.pi)*mp.e**(-3*iu*mp.pi/4)*mp.sqrt(eps)*mp.e**(-mp.pi/(2*eps))
    Iv = I_complex(iu*eps)
    print(f"eps={mp.nstr(eps,6)}:")
    print(f"   |loop_H|/|Jpred| = {mp.nstr(abs(lp)/abs(Jpred), 12)}   arg(loop_H)-arg(pred) = {mp.nstr(mp.arg(lp)-mp.arg(Jpred), 8)}")
    print(f"   |S_+^R - S_-^R| = {mp.nstr(abs(SRp-SRm),4)} (should be ~e^{-pi/eps} = {mp.nstr(mp.e**(-mp.pi/eps),3)})")
    print(f"   |J|={mp.nstr(abs(J),10)} |Jpred|={mp.nstr(abs(Jpred),10)} ratio={mp.nstr(abs(J)/abs(Jpred),10)} dphase={mp.nstr(mp.arg(J)-mp.arg(Jpred),7)}")
    print(f"   |S_- - I(i eps)| = {mp.nstr(abs(Sm-Iv),3)}")
    res[str(eps)] = (mp.nstr(J, 30), mp.nstr(Jpred, 30))
json.dump(res, open("stokes_jump.json","w"))
