# 03_pieces.py  -- split Q0 = 6 (P1 - P2a - P2c) and evaluate each piece to high precision.
#   P1  = int_0^1 [2 log(1+q)/(1-q^2) - log(2q)/(1-q)] ds                     (2-D "squircle" world)
#   P2a = int_0^1 asinh(q)/sqrt(1+q^2) ds,  P2c = int_0^1 asinh(1/q)/sqrt(1+q^2) ds   (vertex world)
# with q = sqrt(s^2+(1-s)^2).  Also the transformed 1-D forms used in the proofs:
#   P2a = (sqrt2/4) int_1^{x1} log(Y) (1-1/Y) dx/x ,  Y = (sqrt3 x^2+4x+sqrt3)/(2x),  x1 = sqrt3 (1+sqrt2)
#   P2c = int_{sqrt(3/2)}^{sqrt2} log((beta+1)/(beta-1)) dbeta/sqrt(2 beta^2-3)
import mpmath as mp, sys, json
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 120
mp.mp.dps = DPS + 20

def with_m(fun):
    def g(s):
        s = mp.mpf(s); m = s*(1-s)
        extra = int(max(0, -mp.log10(m))) + 10 if m > 0 else 0
        with mp.workdps(mp.mp.dps + extra):
            q = mp.sqrt(1-2*m)
            v = fun(q, m)
        return +v
    return g

P1 = 2*mp.quad(with_m(lambda q, m: (2*mp.log(1+q) - (1+q)*mp.log(2*q))/(2*m)), [0, mp.mpf(1)/4, mp.mpf(1)/2])
P2a = 2*mp.quad(with_m(lambda q, m: mp.asinh(q)/mp.sqrt(1+q*q)), [0, mp.mpf(1)/2])
P2c = 2*mp.quad(with_m(lambda q, m: mp.asinh(1/q)/mp.sqrt(1+q*q)), [0, mp.mpf(1)/2])
s3, s2 = mp.sqrt(3), mp.sqrt(2)
x1 = s3*(1+s2)
P2a_t = s2/4*mp.quad(lambda x: mp.log((s3*x*x+4*x+s3)/(2*x))*(1-2*x/(s3*x*x+4*x+s3))/x, [1, 2, x1])
P2c_t = mp.quad(lambda be: mp.log((be+1)/(be-1))/mp.sqrt(2*be*be-3), [mp.sqrt(mp.mpf(3)/2), s2])
Q0 = 6*(P1 - P2a - P2c)
print("P1  =", mp.nstr(P1, DPS))
print("P2a =", mp.nstr(P2a, DPS), "  transformed-form residual", mp.nstr(P2a-P2a_t, 3))
print("P2c =", mp.nstr(P2c, DPS), "  transformed-form residual", mp.nstr(P2c-P2c_t, 3))
print("Q0  = 6(P1-P2a-P2c) =", mp.nstr(Q0, DPS))
try:
    Q0k = mp.mpf(open("Q0_value.txt").read().strip())
    print("   vs k^3 representation: residual", mp.nstr(Q0-Q0k, 3))
except Exception:
    pass
json.dump({"P1": mp.nstr(P1, DPS+5), "P2a": mp.nstr(P2a, DPS+5), "P2c": mp.nstr(P2c, DPS+5),
           "Q0": mp.nstr(Q0, DPS+5)}, open("pieces_%d.json" % DPS, "w"), indent=1)
