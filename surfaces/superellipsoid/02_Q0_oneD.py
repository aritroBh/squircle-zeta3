# 02_Q0_oneD.py
# Checks the 1-D reduction of the vertex constant
#   Q0 = int_Delta [sqrt(u^2+v^2)+sqrt(v^2+w^2)+sqrt(w^2+u^2) - sqrt(u^2+v^2+w^2) - 1]/(uvw) dudv
#      = 6 * int_0^1 J(q(s)) ds,   q = sqrt(s^2+(1-s)^2),
#   J(q) = [2 log(1+q) - (1+q) log(2q)]/(1-q^2) - [asinh(q)+asinh(1/q)]/sqrt(1+q^2)
# against (a) the k^3 representation (01_kfun_Q0.py, Q0_value.txt) and (b) a direct 2-D
# quadrature of the simplex integral (lower precision, independent).
import mpmath as mp, sys
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.mp.dps = DPS + 15

def Js(s):
    """J(q(s)) evaluated with extra working precision near s=0 (removable 0/0)."""
    s = mp.mpf(s)
    m = s*(1-s)
    extra = int(max(0, -mp.log10(m))) + 10 if m > 0 else 0
    with mp.workdps(mp.mp.dps + extra):
        q = mp.sqrt(1 - 2*m)
        first = (2*mp.log(1+q) - (1+q)*mp.log(2*q))/(2*m)
        val = first - (mp.asinh(q) + mp.asinh(1/q))/mp.sqrt(1 + q*q)
    return +val

Q0_1d = 12*mp.quad(Js, [0, mp.mpf(1)/4, mp.mpf(1)/2])
print("Q0 (1-D J-formula)      =", mp.nstr(Q0_1d, DPS))
try:
    Q0_k = mp.mpf(open("Q0_value.txt").read().strip())
    print("Q0 (k^3 representation) =", mp.nstr(Q0_k, DPS))
    print("residual                =", mp.nstr(Q0_1d - Q0_k, 3))
except Exception as e:
    print("no Q0_value.txt", e)

# (b) direct 2-D quadrature of R_0 over the simplex (fundamental region u>=v>=w, x6),
# polar coordinates about the corner (1,0,0): v = r t, w = r (1-t), t in [1/2,1], r in [0,1/(1+t)]
if "--direct" in sys.argv:
    mp.mp.dps = 25
    def R0(u, v, w):
        return (mp.sqrt(u*u+v*v) + mp.sqrt(v*v+w*w) + mp.sqrt(w*w+u*u)
                - mp.sqrt(u*u+v*v+w*w) - 1)/(u*v*w)
    def inner(t):
        return mp.quad(lambda r: r*R0(1-r, r*t, r*(1-t)), [0, 1/(1+t)])
    Qd = 6*mp.quad(inner, [mp.mpf(1)/2, 1])
    print("Q0 (direct 2-D simplex) =", mp.nstr(Qd, 20), " diff", mp.nstr(Qd - Q0_1d, 3))
