# SM6 item 1: c6, c7 from the SM5 pickled quadrature values.
# Y(tau) = (I - 8 - pi^2/3 t^2 - 2z3 t^3 - z4/4 t^4 - c5 t^5)/t^6 = c6 + c7 t + c8 t^2 + ...
import mpmath as mp, pickle, sys
mp.mp.dps = 170
data = pickle.load(open("../SM5_softmax_c5/sm5_I_values.pkl", "rb"))
z = mp.zeta
c5 = mp.mpf(107)/32*z(5) - z(2)*z(3)
def run(i0, label):
    ds = data[i0:]
    taus = [mp.mpf(1)/n for n, _, _ in ds]
    I = [mp.mpf(v) for _, v, _ in ds]
    Y = [(I[i]-8-mp.pi**2/3*taus[i]**2-2*z(3)*taus[i]**3-z(4)/4*taus[i]**4-c5*taus[i]**5)/taus[i]**6 for i in range(len(I))]
    T = [Y[:]]
    for k in range(1, len(Y)):
        f = mp.mpf(2)**k; p = T[-1]
        T.append([(f*p[i+1]-p[i])/(f-1) for i in range(len(p)-1)])
    # c6 = T[k][-1]; for c7, apply Richardson to (Y-c6)/tau is noisy; instead eliminate tau^1 with weights: second sequence
    # c7 estimate: Z = (Y - c6_est)/tau at each level using best c6, then Neville again
    print("== %s: points n=%d..%d (%d pts)" % (label, ds[0][0], ds[-1][0], len(ds)))
    for k in range(4, len(T)):
        d = abs(T[k][-1]-T[k-1][-1]) if k>0 else 0
        print("k=%2d c6 = %s   |diff prev order| = %s" % (k, mp.nstr(T[k][-1], 50), mp.nstr(d, 3)))
    return T, taus, Y
if __name__ == "__main__":
    for i0, lab in [(0, "all"), (3, "n>=160"), (5, "n>=640")]:
        run(i0, lab)
