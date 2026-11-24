# 14_Q1_pieces.py -- the first-order vertex coefficient Q1 = Q'(0) from the representation
#   Q(eps) = -(1-eps)^(-2) int_Delta D(sigma) prod sigma_i^(theta-1) (sum sigma_i^(1+theta))^(-1-2eps) dsigma,
#   theta = eps/(1-eps),  D(x) = sum x_i - sum_{i<j} |x_ij| + |x|,
# which gives  Q1 = 2 Q0 - int_Delta D/(s1 s2 s3) * sum_i (1-s_i) log s_i
#                 = 2 Q0 - 3 Zq + 6 Za + 6 Zb, with
#   Psi1(c) = int_0^1 [sqrt(u^2 + c^2 (1-u)^2) - c (1-u)] log(u) du/u,
#   Zq = int_0^1 [Psi1(q(s)) - Psi1(1)] ds/(s(1-s)),  Za = int_0^1 [Psi1(t)-Psi1(1)] dt/(1-t),
#   Zb = int_0^1 [Psi1(t)+1] dt/t.
import mpmath as mp, sys, json
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
mp.mp.dps = DPS + 20

def Psi1(c):
    c = mp.mpf(c)
    f = lambda u: (mp.sqrt(u*u + c*c*(1-u)**2) - c*(1-u))*mp.log(u)/u
    # split near the scale where u ~ c
    pts = [0, min(c, mp.mpf(1)/2)/4, min(c, mp.mpf(1)/2), 1] if c > 0 else [0, 1]
    pts = sorted(set(pts))
    return mp.quad(f, pts)

if __name__ == "__main__":
    P1v = Psi1(1)
    print("Psi1(0) =", mp.nstr(Psi1(0), 20), "(expect -1)")
    def Zq_int(s):
        m = s*(1-s)
        with mp.workdps(mp.mp.dps + int(max(0, -mp.log10(m))) + 5):
            q = mp.sqrt(1 - 2*m)
            v = (Psi1(q) - P1v)/m
        return +v
    Zq = 2*mp.quad(Zq_int, [0, mp.mpf(1)/8, mp.mpf(1)/2])
    Za = mp.quad(lambda t: (Psi1(t) - P1v)/(1-t), [0, mp.mpf(1)/8, mp.mpf(1)/2, 1])
    Zb = mp.quad(lambda t: (Psi1(t) + 1)/t, [0, mp.mpf(1)/64, mp.mpf(1)/8, mp.mpf(1)/2, 1])
    Q0 = mp.mpf(open("Q0_value.txt").read().strip())
    Q1 = 2*Q0 - 3*Zq + 6*Za + 6*Zb
    print("Psi1(1) =", mp.nstr(P1v, DPS))
    print("Zq =", mp.nstr(Zq, DPS)); print("Za =", mp.nstr(Za, DPS)); print("Zb =", mp.nstr(Zb, DPS))
    print("Q1 =", mp.nstr(Q1, DPS))
    for fn in ("Q1_value_105.txt", "Q1_value_50.txt"):
        try:
            ref = mp.mpf(open(fn).read().strip()); print("  vs %s: %s" % (fn, mp.nstr(Q1-ref, 3))); break
        except Exception: pass
    json.dump({k: mp.nstr(v, DPS) for k, v in dict(Psi1_1=P1v, Zq=Zq, Za=Za, Zb=Zb, Q1=Q1).items()},
              open("Q1_pieces_%d.json" % DPS, "w"), indent=1)
