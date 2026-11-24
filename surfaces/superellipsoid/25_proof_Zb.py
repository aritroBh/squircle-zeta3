# 25_proof_Zb.py -- proof of  Zb = 2 M1 - 2 M0 + pi^2/12 - log^2(2)/2.
# Reduction (see result.md):  Zb = -int_{x+y<=1} D2(x,y) log y dx dy/(xy) = F(0) - X1 - X2,
#   F(0) = -2 M0,  X2 = int_0^1 D2(s,1-s) log s ds/s = (F'(0) - F(0))/2 = M0 - M1,
#   X1 = int_0^1 (1-q) log(1-s) ds/s,  q = sqrt(s^2+(1-s)^2).
# So the claim is equivalent to  X1 = -M1 - M0 - pi^2/12 + log^2(2)/2,  proved here exactly:
# Euler substitution s=(v-b)(v+a)/(4v):  X1 = (1/(2 sqrt2)) int_b^a R(v) [log(a-v)+log(v+b)-log v-2 log2] dv,
#   R = 1/v^2 - 2b/v + 4 sqrt2/(v+a) - 1,  integrated by parts against P(v) = int_a^v R  (P(a)=0).
import mpmath as mp
from fractions import Fraction as Fr
from qs2 import Q2, W2, Li2, PI2, R2, A, B, logprod, logsq, logw, ONEW, canon_zero_li2
from dilog_engine import pool, Reducer
mp.mp.dps = 50
r2 = mp.sqrt(2); a = 1+r2; b = r2-1; lam = mp.log(2); L = mp.log(a)
def chk(name, w, val):
    d = w.num() - val
    print("  %-40s check %s" % (name, mp.nstr(d, 3))); assert abs(d) < mp.mpf(10)**-40, name

def I(p, r, Aa, Bb):   # int_A^B log(v+p) dv/(v+r), v+p>0, v+r>0 on [A,B]
    p, r, Aa, Bb = Q2(1)*p, Q2(1)*r, Q2(1)*Aa, Q2(1)*Bb
    if p == r: return logsq(Bb+p, Fr(1, 2)) + logsq(Aa+p, Fr(-1, 2))
    c = p - r; ylo, yhi = Aa + r, Bb + r
    if float(c) > 0: return logprod(c, yhi) - logprod(c, ylo) - Li2(-yhi/c) + Li2(-ylo/c)
    d = -c
    return logsq(yhi, Fr(1, 2)) - logsq(ylo, Fr(1, 2)) + Li2(d/yhi) - Li2(d/ylo)
def Ireg(p, c, Aa, Bb):  # int_A^B log((v+p)/(c+p)) dv/(v-c), allowed with c an endpoint
    p, c, Aa, Bb = Q2(1)*p, Q2(1)*c, Q2(1)*Aa, Q2(1)*Bb
    return Li2(-(Aa-c)/(c+p)) - Li2(-(Bb-c)/(c+p))

lo, hi = B, A
Pb = logw(Q2(2), 2*R2) + logw(A, -4)               # P(b) = 2 sqrt2 log2 - 4L
Pfun = lambda v: -1/v + 1/a - v + a - 2*b*mp.log(v/a) + 4*r2*mp.log((v+a)/(2*a))
chk("P(b)", Pb, Pfun(b))
# int_b^a Prat/(v-c) for c = a, -b, 0  (Prat = -1/v + 1/a - v + a)
PratInt = {"a": logw(A, 2*B) + ONEW(-2), "-b": logw(Q2(2), 2*R2) + logw(A, 2*B) + ONEW(-2), "0": logw(A, 4*R2) + ONEW(-4)}
for key, cc in (("a", a), ("-b", -b), ("0", 0)):
    chk("Prat/(v-c), c=%s" % key, PratInt[key], mp.quad(lambda v: (-1/v + 1/a - v + a)/(v - cc), [b, a]))
# log parts of P:  -2b log(v/a) + 4 sqrt2 log((v+a)/(2a))
def Plog_int(c):
    """ int_b^a [ -2b log(v/a) + 4sqrt2 log((v+a)/(2a)) ] dv/(v-c) """
    if c == "a":
        return Ireg(0, A, lo, hi).scale(-2*B) + Ireg(A, A, lo, hi).scale(4*R2)
    if c == "-b":   # 1/(v+b): use I with r=b
        t1 = I(0, B, lo, hi) - (logprod(A, hi + B) - logprod(A, lo + B))            # int log(v/a)/(v+b)
        t2 = I(A, B, lo, hi) - (logprod(2*A, hi + B) - logprod(2*A, lo + B))        # int log((v+a)/(2a))/(v+b)
        return t1.scale(-2*B) + t2.scale(4*R2)
    if c == "0":
        t1 = logsq(hi/A, Fr(1, 2)) - logsq(lo/A, Fr(1, 2))
        t2 = I(A, 0, lo, hi) - (logprod(2*A, hi) - logprod(2*A, lo))
        return t1.scale(-2*B) + t2.scale(4*R2)
for key, cc in (("a", a), ("-b", -b), ("0", 0)):
    chk("Plog/(v-c), c=%s" % key, Plog_int(key), mp.quad(lambda v: (-2*b*mp.log(v/a) + 4*r2*mp.log((v+a)/(2*a)))/(v - cc), [b, a]))
PI_ = {k: PratInt[k] + Plog_int(k) for k in ("a", "-b", "0")}
# int_b^a R*l dv = -P(b) l(b) - int P l'   for l = log(a-v) [l' = 1/(v-a)], log(v+b) [1/(v+b)], log v [1/v]; const: -P(b)*const
T_amv = logprod(Q2(2), Q2(2), 0)  # placeholder zero
def pl(x):  # -P(b)*log x as W2:  P(b) is weight 1 -> product
    # P(b) = 2sqrt2 log2 - 4L  ->  -P(b) log x = -2sqrt2 log2 log x + 4 L log x
    return logprod(Q2(2), x, -2*R2) + logprod(A, x, 4)
T1 = pl(A - B) - PI_["a"]          # log(a-v), l(b) = log(a-b) = log 2
T2 = pl(2*B) - PI_["-b"]           # log(v+b), l(b) = log 2b
T3 = pl(B) - PI_["0"]              # log v,    l(b) = log b
T4 = (logw(Q2(2), 2*R2) + logw(A, -4))  # P(b) ; const term = -2 log2 * (-P(b))
T4 = logprod(Q2(2), Q2(2), 4*R2) + logprod(Q2(2), A, -8)   # = -(-2log2)*P(b) ... see below
# const contribution: int R * (-2 log2) = -2log2 * (P(a)-P(b)) = 2 log2 * P(b)
Rf = lambda v: 1/v**2 - 2*b/v + 4*r2/(v+a) - 1
chk("int R log(a-v)", T1, mp.quad(lambda v: Rf(v)*mp.log(a-v), [b, a]))
chk("int R log(v+b)", T2, mp.quad(lambda v: Rf(v)*mp.log(v+b), [b, a]))
chk("int R log v", T3, mp.quad(lambda v: Rf(v)*mp.log(v), [b, a]))
chk("int R (-2log2)", T4, mp.quad(lambda v: Rf(v)*(-2*lam), [b, a]))
X1 = (T1 + T2 - T3 + T4).scale(Q2(0, Fr(1, 4)))      # 1/(2 sqrt2) = sqrt2/4
X1num = mp.quad(lambda s: (1 - mp.sqrt(s*s+(1-s)**2))*mp.log(1-s)/s, [0, mp.mpf(1)/2, 1])
chk("X1 (exact, unreduced)", X1, X1num)
red = Reducer(pool(), basis=[B, -(B**4), Q2(Fr(1, 3))])
X1r, cert = red.reduce(canon_zero_li2(X1))
print("X1 =", X1r)
# target: -M1 - M0 - pi^2/12 + log^2 2/2, with M0 = log2 - sqrt2 L,
# M1 = (19sqrt2/48 - 1/24) pi^2 + log^2 2/2 + 2 sqrt2 log2 L - (1+5sqrt2/2) L^2 - 5 sqrt2 Li2(b)
M0w = logw(Q2(2)) + logw(A, -R2)
M1w = (PI2(Q2(Fr(-1, 24), Fr(19, 48))) + logsq(Q2(2), Fr(1, 2)) + logprod(Q2(2), A, 2*R2)
       + logsq(A, Q2(-1, Fr(-5, 2))) + Li2(B, Q2(0, -5)))
target = -M1w - M0w + PI2(Fr(-1, 12)) + logsq(Q2(2), Fr(1, 2))
diff, _ = red.reduce(canon_zero_li2(X1r - target))
diff = canon_zero_li2(diff)
print("X1 - (-M1 - M0 - pi^2/12 + log^2 2/2)  (exact) =", diff if diff.d else "0")
if not diff.d:
    print("hence Zb = M1 - 3 M0 - X1 = 2M1 - 2M0 + pi^2/12 - log^2(2)/2   (PROVED)")
