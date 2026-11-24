# 17_Q1_subpieces.py -- elementary 1-D integrals for Q1:
#   Q1 = 2 Q0 - 3 Zq + 6 Za + 6 Zb,   Zq = Zq1 + Zq2 + Zq3,
#   Zq1 = -sqrt2 [ 2 V3' - ... ]  (see below),  with q_u = sqrt(u^2+(1-u)^2):
#   Phi(u) = -sqrt2 (1-u) log X - 2 q_u log Y - 2 (1-u)(log2 - sqrt2 L),
#   X = (sqrt2 q_u + 1-u)/(sqrt2 q_u - 1+u),  Y = (3u^2-2u+1)/(2 q_u^2),   Zq = int_0^1 Phi(u) log u du/u.
#   V1 = int q_u log(3u^2-2u+1) log u du/u
#   V2 = int (1-u) log(3u^2-2u+1) log u du/u        (closed form below)
#   V3 = int [2(1-u) log(sqrt2 q_u+1-u) - 2L] log u du/u
#   V4 = int [q_u log(2 q_u^2) - log2] log u du/u
#   Zq = -sqrt2 (V3 - V2) - 2 V1 + 2 V4 - 2(log2 - sqrt2 L)
#   Zb = int_0^1 [q_u - 1 - u log((u+q_u)/(2u))] log u du/u     (swap of integrations)
import mpmath as mp, sys, json
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.mp.dps = DPS + 25
r2 = mp.sqrt(2); L = mp.log(1+r2); lam = mp.log(2)
def qu(u): return mp.sqrt(u*u + (1-u)**2)
def I(f):  # integrals over (0,1) with log singularities at 0
    return mp.quad(f, [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1])
def safe(f):
    def g(u):
        u = mp.mpf(u)
        with mp.workdps(mp.mp.dps + (int(-mp.log10(u)) + 10 if u < mp.mpf('0.01') else 0)):
            v = f(u)
        return +v
    return g
V1 = I(safe(lambda u: qu(u)*mp.log(3*u*u-2*u+1)*mp.log(u)/u))
V2 = I(safe(lambda u: (1-u)*mp.log(3*u*u-2*u+1)*mp.log(u)/u))
V3 = I(safe(lambda u: (2*(1-u)*mp.log(r2*qu(u)+1-u) - 2*L)*mp.log(u)/u))
V4 = I(safe(lambda u: (qu(u)*mp.log(2*qu(u)**2) - lam)*mp.log(u)/u))
Zq = -r2*(V3 - V2) - 2*V1 + 2*V4 - 2*(lam - r2*L)
Zb = I(safe(lambda u: (qu(u) - 1 - u*mp.log((u+qu(u))/(2*u)))*mp.log(u)/u))
#   Za = int_0^1 [1 - q + (1-u) log((q+u-1)/u) + q log(2q^2/(u(u+q)))] log u du/u
def za_int(u):
    q = qu(u)
    # q+u-1 = (q^2-(1-u)^2)/(q+1-u) = u^2/(q+1-u)  (avoids cancellation)
    return (1 - q + (1-u)*mp.log(u/(q+1-u)) + q*mp.log(2*q*q/(u*(u+q))))*mp.log(u)/u
Za = I(safe(za_int))
w = 1 - 1j*r2
V2c = 2*mp.re(mp.polylog(3, w)) - 2*mp.re(-mp.log(1-w) + mp.log(1-w)/w + 2 - mp.polylog(2, w)/w)
print("V2 closed form residual:", mp.nstr(V2 - V2c, 3))
Q0 = mp.mpf(json.load(open('Q1_hp_60.json'))['Q0']) if DPS<=60 else None
for k, v in dict(V1=V1, V2=V2, V3=V3, V4=V4, Zq=Zq, Zb=Zb, Za=Za).items():
    print(k, "=", mp.nstr(v, DPS))
try:
    ref = json.load(open("Q1_hp_60.json"))
    print("Zq vs Psi1-route:", mp.nstr(Zq - mp.mpf(ref["Zq"]), 3), "  Zb:", mp.nstr(Zb - mp.mpf(ref["Zb"]), 3), "  Za:", mp.nstr(Za - mp.mpf(ref["Za"]), 3))
except Exception as e:
    print(e)
json.dump({k: mp.nstr(v, DPS) for k, v in dict(V1=V1, V2=V2, V3=V3, V4=V4, Zq=Zq, Zb=Zb, Za=Za).items()},
          open("Q1_sub_%d.json" % DPS, "w"), indent=1)
