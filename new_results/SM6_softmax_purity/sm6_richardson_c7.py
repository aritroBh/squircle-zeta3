# SM6 item 1 (c7): Z(tau) = (Y - c6)/tau = c7 + c8 tau + ..., Neville ratio 2. Two choices of c6.
import mpmath as mp, pickle
mp.mp.dps = 170
data = pickle.load(open("../SM5_softmax_c5/sm5_I_values.pkl", "rb"))
z = mp.zeta
c5 = mp.mpf(107)/32*z(5) - z(2)*z(3)
taus = [mp.mpf(1)/n for n, _, _ in data]
I = [mp.mpf(v) for _, v, _ in data]
Y = [(I[i]-8-mp.pi**2/3*taus[i]**2-2*z(3)*taus[i]**3-z(4)/4*taus[i]**4-c5*taus[i]**5)/taus[i]**6 for i in range(len(I))]
def neville(Y):
    T = [Y[:]]
    for k in range(1, len(Y)):
        f = mp.mpf(2)**k; p = T[-1]
        T.append([(f*p[i+1]-p[i])/(f-1) for i in range(len(p)-1)])
    return T
T = neville(Y)
c6_rich = T[-1][0]
print("c6 (Richardson, best) =", mp.nstr(c6_rich, 60))
def show(c6, label):
    Zs = [(Y[i]-c6)/taus[i] for i in range(len(Y))]
    T2 = neville(Zs)
    print("== c7 using", label)
    for k in range(4, len(T2)):
        print("k=%2d c7 = %s |diff prev| = %s" % (k, mp.nstr(T2[k][-1], 50), mp.nstr(abs(T2[k][-1]-T2[k-1][-1]), 3)))
    return T2[-1][0]
c7a = show(c6_rich, "Richardson c6")
mp.mp.dps = 170
# exact-gamma c6 (quadrature at 80 digits)
c6g = mp.mpf("0.042950827880284033929433804835135265510296727692058")
c7b = show(c6g, "c6 from gamma formula (52 digits)")
open("sm6_c7_numeric.txt", "w").write(mp.nstr(c7a, 60)+"\n")
