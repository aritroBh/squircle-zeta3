# SM6 item 4 numerics: p-norm perimeter from the exact 1-D integral, Richardson for c6, c7, compare with the pure closed forms.
import mpmath as mp
from multiprocessing import Pool
DPS = 150
def I_of(n):
    mp.mp.dps = DPS
    tau = mp.mpf(1)/n
    w = lambda v: mp.sqrt(1+mp.tanh(2*v))/mp.cosh(2*v)
    f = lambda v: w(v)*((1+mp.exp(-2*v/tau))**tau-1)
    pts = [mp.mpf(0)]; k = mp.mpf(1)
    top = 160*tau
    while k*tau < top: pts.append(k*tau); k *= 2
    pts.append(top)
    val = mp.quad(f, pts, maxdegree=9)
    return n, 4*mp.sqrt(2)*(1+val)
if __name__ == "__main__":
    ns = [10*2**i for i in range(15)]
    with Pool(8) as P: res = P.map(I_of, ns)
    mp.mp.dps = DPS
    z = mp.zeta; s2 = mp.sqrt(2)
    c = {0: 4*s2, 2: s2*z(2), 3: s2*z(3), 4: -33*s2/16*z(4), 5: -s2*(mp.mpf(3)/4*z(5)+mp.mpf(7)/4*z(2)*z(3)),
         6: s2*(mp.mpf(4039)/256*z(6)-mp.mpf(7)/8*z(3)**2),
         7: s2*(mp.mpf(1071)/64*z(4)*z(3)+mp.mpf(141)/16*z(2)*z(5)+mp.mpf(129)/16*z(7))}
    taus = [mp.mpf(1)/n for n, _ in res]; I = [v for _, v in res]
    def neville(Y):
        T = [Y[:]]
        for k in range(1, len(Y)):
            f = mp.mpf(2)**k; p = T[-1]
            T.append([(f*p[i+1]-p[i])/(f-1) for i in range(len(p)-1)])
        return T
    Y6 = [(I[i]-sum(c[m]*taus[i]**m for m in range(0, 6) if m in c))/taus[i]**6 for i in range(len(I))]
    T = neville(Y6)
    print("squircle c6 formula:", mp.nstr(c[6], 40))
    for k in (6, 8, 10, 12, 14):
        print(" k=%2d Richardson c6 = %s  diff to formula %s" % (k, mp.nstr(T[k][-1], 40), mp.nstr(T[k][-1]-c[6], 3)))
    Z7 = [(Y6[i]-c[6])/taus[i] for i in range(len(I))]
    T2 = neville(Z7)
    print("squircle c7 formula:", mp.nstr(c[7], 40))
    for k in (6, 8, 10, 12, 14):
        print(" k=%2d Richardson c7 = %s  diff to formula %s" % (k, mp.nstr(T2[k][-1], 40), mp.nstr(T2[k][-1]-c[7], 3)))
