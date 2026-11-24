# SM5 item 1: c5 of I(tau) = int_0^{2pi} h dtheta (soft-max of the square), dps 110.
import mpmath as mp, sys, pickle
from multiprocessing import Pool
DPS = 170
def pts(b, tau):
    p = [mp.mpf(0)]; k = mp.mpf(1)
    while k*tau < b: p.append(k*tau); k *= 2
    p.append(b)
    return p                       # breakpoints geometric from theta=0 (layer at theta_j=0 side)
def I_of(n):
    mp.mp.dps = DPS
    tau = mp.mpf(1)/n
    # use symmetry: I = 8 * int_0^{pi/4}; layer sits at the kink theta=0 (and by symmetry pi/4 is not a kink)
    def f(t):
        c, s = mp.cos(t), mp.sin(t)
        return (c+s)*(1+mp.exp(-2*s/tau)+mp.exp(-2*c/tau)+mp.exp(-2*(c+s)/tau))**tau
    v, e = mp.quad(f, pts(mp.pi/4, tau), error=True, maxdegree=12)
    return n, 8*v, 8*e
if __name__ == "__main__":
    ns = [20*2**k for k in range(18)]
    with Pool(14) as P: res = P.map(I_of, ns)
    mp.mp.dps = DPS; pickle.dump([(n, mp.nstr(v, DPS-5), mp.nstr(e, 5)) for n, v, e in res], open("sm5_I_values.pkl", "wb"))
    for n, v, e in res: print(n, mp.nstr(e, 3))
