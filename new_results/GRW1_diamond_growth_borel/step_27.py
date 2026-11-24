
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
iu = mp.mpc(0,1)
alpha = mp.pi/6

def ray_int(sign, eps, fun):
    R = 30*eps
    z_of_r = lambda r: r*mp.e**(iu*(mp.pi/2+sign*alpha))
    w = lambda r: mp.e**(-(r/eps)*mp.e**(sign*iu*alpha))
    pts = [0, eps/2, 2*eps, 6*eps, 12*eps, 20*eps, R]
    return mp.quad(lambda r: w(r)*fun(z_of_r(r)), pts)

def lateral(sign, eps, fun):
    return mp.e**(sign*iu*alpha)/eps * ray_int(sign, eps, fun)

def loop_H_over_eta(eps):
    f = lambda s: mp.e**(-(mp.pi/2+s)/eps)*mp.sqrt(mp.sin(s))*mp.e**(iu*(mp.pi/4 - s/2))
    v = mp.quad(f, [0, eps, 4*eps, mp.pi/2, mp.pi])
    return -4*iu*v/(iu*eps)

def I_complex(eta):
    f = lambda u: mp.exp(-u)*Hfun(eta*u)*(1+mp.exp(-u))**(-2-eta)
    return mp.quad(f, [0, 1, 4, 12, 40, 120])

print("Stokes jump: full J = loop_H/(i eps) + (S_+^R - S_-^R) vs pred")
res = {}
for eps in [mp.mpf('0.2'), mp.mpf('0.15'), mp.mpf('0.1'), mp.mpf('0.075')]:
    SRp = lateral(+1, eps, pevalR); SRm = lateral(-1, eps, pevalR)
    SHm = lateral(-1, eps, Hfun)
    lp = loop_H_over_eta(eps)
    Sm = SHm + SRm
    J = lp + (SRp - SRm)
    Jpred = 2*mp.sqrt(mp.pi)*mp.e**(-3*iu*mp.pi/4)*mp.sqrt(eps)*mp.e**(-mp.pi/(2*eps))
    Iv = I_complex(iu*eps)
    print(f"eps={mp.nstr(eps,6)}:")
    print(f"   |loop|/|pred|={mp.nstr(abs(lp)/abs(Jpred),10)}  dphase_loop={mp.nstr(mp.arg(lp)-mp.arg(Jpred),8)}")
    print(f"   |S_+^R-S_-^R|={mp.nstr(abs(SRp-SRm),3)}   |J|/|pred|={mp.nstr(abs(J)/abs(Jpred),10)}  dphase_J={mp.nstr(mp.arg(J)-mp.arg(Jpred),8)}")
    print(f"   |S_- - I(i eps)|={mp.nstr(abs(Sm-Iv),3)}")
    res[str(eps)] = (mp.nstr(J, 30), mp.nstr(Jpred, 30))
json.dump(res, open("stokes_jump.json","w"))
