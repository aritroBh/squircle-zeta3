import pickle, sys, time
sys.path.insert(0, "<project>/research/08_theory/E5c_nonlocal_z5")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 60
import symeng
r3 = mp.sqrt(3)
M = []
for j in range(6):
    c, s = mp.cos(j*mp.pi/6), mp.sin(j*mp.pi/6)
    M.append(((1 if j<=3 else -1)*c, (1 if j<=3 else -1)*s))
y0 = (2-r3, mp.mpf(0)); y1 = ((2*r3-3)/2, (2-r3)/2)
P = symeng.edge_params(M, y0, y1)
t0 = time.time()
I = mp.quad(lambda t: symeng.E5_sym(t, M, P, conv='B'), [0, 1])
P5_engine = mp.mpf('-62.5372477549863953530912928220580357268476927741037326268628')
print("1440 * int E5 =", mp.nstr(1440*I, 40))
print("P5 engine     =", mp.nstr(P5_engine, 40))
print("diff:", mp.nstr(abs(1440*I - P5_engine), 4), " (%.1fs)" % (time.time()-t0))
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
mons = {}
for (n, le), c in E5X.items():
    key = (n, le[0], le[1], le[2], le[3], le[4], le[5])
    constdeg = le[6]+le[7]+le[8]
    mons.setdefault(key, 0)
print("distinct (n,LT,LU,W0,W1,W2,W3) monomials:", len(mons))
w5 = [k for k in mons if sum(k[1:]) <= 5]
print("with variable weight <=5:", len(w5))
from collections import Counter
print("by variable weight:", dict(sorted(Counter(sum(k[1:]) for k in mons).items())))
