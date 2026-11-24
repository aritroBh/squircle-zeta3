# Planner check of K1e: P'''(1) of the fresh hexagon V0=(1,0), V1=(3,3), V2=(0,1) (+negatives), weights w=(3,5,3),
# computed directly from the polar perimeter integral of K_p = {sum_j w_j |k_j cross(V_j,x)|^p <= 1} (k_j = 1/18 from N(V_q)=1),
# with no use of K1e's engine. Central finite differences in p at p = 1, Richardson-extrapolated.
import mpmath as mp, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 80
V = [(1, 0), (3, 3), (0, 1)]; w = [3, 5, 3]; k = [mp.mpf(1)/18]*3
# k from N(V_q) = 1: check
for q in range(3):
    assert abs(sum(w[j]*k[j]*abs(V[j][0]*V[q][1]-V[j][1]*V[q][0]) for j in range(3)) - 1) < mp.mpf(10)**(-50)
def P(p):
    c = [w[j]*k[j]**p for j in range(3)]
    def integrand(th):
        s, co = mp.sin(th), mp.cos(th)
        f = 0; fp = 0
        for j in range(3):
            sj = V[j][0]*s - V[j][1]*co; dsj = V[j][0]*co + V[j][1]*s
            a = abs(sj)
            f += c[j]*a**p; fp += c[j]*p*a**(p-1)*mp.sign(sj)*dsj
        rho = f**(-1/p); drho = -(1/p)*f**(-1/p-1)*fp
        return mp.sqrt(rho**2 + drho**2)
    br = [0, mp.pi/4, mp.pi/2, mp.pi]
    return 2*sum(mp.quad(integrand, [br[i], br[i+1]], method='tanh-sinh') for i in range(3))
print("P(1) =", mp.nstr(P(mp.mpf(1)), 30))
def d3(h):
    return (P(1+2*h) - 2*P(1+h) + 2*P(1-h) - P(1-2*h))/(2*h**3)
h = mp.mpf(10)**(-6)
a, b = d3(h), d3(h/2)
rich = (4*b - a)/3
print("d3(h) =", mp.nstr(a, 30)); print("d3(h/2) =", mp.nstr(b, 30)); print("Richardson P'''(1) =", mp.nstr(rich, 30))
ref = mp.mpf("-44.5799197716825334100277211047")
print("K1e formula value  =", ref, " diff =", mp.nstr(rich-ref, 3))
