"""verify step: the definition of M_k against the actual squircle perimeter.
P(p)/8 = (eps/2) int_0^1 sqrt(s^(2eps-2)+(1-s)^(2eps-2)) ds with x=s^eps, y=(1-s)^eps (eps=1/p), reduced to a smooth integral and computed by mp.quad.
Compare with 1 + sum_{k=0}^{K} M_k / p^(k+1) using M_k from verify_lib.direct_M (own definition-based quadrature).
The remainder must scale like M_{K+1}/p^(K+2).  usage: verify_perimeter.py [K]"""
import sys, os, math
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, W)
import verify_lib as VL
import mpmath as mp
K = int(sys.argv[1]) if len(sys.argv) > 1 else 7
TGT = 30
Y = math.log(10) * (TGT + 12) + 25 + 3 * K
Ms = []
for k in range(K + 2):
    VL.set_dps(TGT + int(Y / 2.302585) + 30)
    Ms.append(VL.direct_M(k, TGT, Y=Y))
VL.set_dps(50)
print("mp.mp.dps =", mp.mp.dps)
print("M_k (own):", [mp.nstr(m, 12) for m in Ms])
for p in (20, 40, 80):
    pp = mp.mpf(p); eps = 1 / pp
    # P/8 = (eps/2) int_0^1 sqrt(s^(2eps-2) + (1-s)^(2eps-2)) ds = int_0^{wmax} sqrt(1 + (s/(1-s))^(2-2eps)) dw,
    # s = w^(1/eps), wmax = (1/2)^eps  (symmetry s <-> 1-s, then s^(eps-1) ds = dw/eps removes the endpoint singularity)
    wmax = mp.mpf(1) / 2 ** eps
    def f(w):
        s_ = w ** pp
        return mp.sqrt(1 + (s_ / (1 - s_)) ** (2 - 2 * eps))
    P8 = mp.quad(f, mp.linspace(0, wmax, 9))
    ser = 1 + sum(Ms[k] / pp ** (k + 1) for k in range(K + 1))
    rem = P8 - ser
    print(f"p={p}: P/8 = {mp.nstr(P8, 25)}; 1+sum_(k<={K}) M_k/p^(k+1) = {mp.nstr(ser, 25)}; remainder = {mp.nstr(rem, 6)}; "
          f"M_(K+1)/p^(K+2) = {mp.nstr(Ms[K+1] / pp ** (K + 2), 6)}  ratio = {mp.nstr(rem / (Ms[K+1] / pp ** (K + 2)), 6)}")
