# SM5 item 1: Richardson extrapolation of R(tau)/tau^5 = c5 + c6 tau + c7 tau^2 + ...
import mpmath as mp, pickle
mp.mp.dps = 170
data = pickle.load(open("sm5_I_values.pkl", "rb"))
z = mp.zeta
taus = [mp.mpf(1)/n for n, _, _ in data]
I = [mp.mpf(v) for _, v, _ in data]
Y = [(I[i] - 8 - mp.pi**2/3*taus[i]**2 - 2*z(3)*taus[i]**3 - z(4)/4*taus[i]**4)/taus[i]**5 for i in range(len(I))]
# Neville in tau (ratio 2): T[k][i] eliminates tau^1..tau^k
T = [Y[:]]
for k in range(1, len(Y)):
    prev = T[-1]; f = mp.mpf(2)**k
    T.append([(f*prev[i+1] - prev[i])/(f-1) for i in range(len(prev)-1)])
print("tau = 1/(20*2^i); column k = k-th Richardson order; entries printed to 40 digits")
out = []
for k, col in enumerate(T):
    s = f"order {k:2d}: " + " | ".join(mp.nstr(c, 40) for c in col[-2:])
    print(s)
best = T[-1][0]
print()
print("diff between last two orders at best level:", mp.nstr(abs(T[-1][0]-T[-2][-1]), 5))
for k in range(5, len(T)):
    print("k=%2d  last entry %s   |last - previous-order last| = %s" % (k, mp.nstr(T[k][-1], 45), mp.nstr(abs(T[k][-1]-T[k-1][-1]), 3)))
mp.mp.dps = 60
open("sm5_c5_numeric.txt", "w").write(mp.nstr(T[-4][-1], 60) + "\n")
