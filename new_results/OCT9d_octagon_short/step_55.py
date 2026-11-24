import sys, time, pickle
import mpmath as mp
sys.path.insert(0, ".")
from oct9d_common import *
from oct9c_common import LAMKEY
t0=time.time()
# the 16 (a,c) pairs from LAMKEY
pairs = {}
for (e, a, c), lbl in LAMKEY.items():
    pairs[lbl] = (sp.Rational(a), sp.Rational(c))
DPS = 460
mp.mp.dps = DPS
lam_hi = {}
for lbl, (a, c) in sorted(pairs.items()):
    am, bm = mp.mpf(str(a)), mp.mpf(str(a*c))
    F = mp.quad(lambda tt: mp.log(tt)*mp.log(1-bm*tt)/(1-am*tt), [0,1])
    lam_hi[lbl] = mp.nstr(am*F, DPS)
    print(lbl, (str(a), str(c)), mp.nstr(am*F, 25), f"({time.time()-t0:.1f}s)")
# consistency vs lam_mine (360)
lam_old = pickle.load(open(base+"OCT9c_octagon_close/lam_mine.pkl","rb"))
for lbl in sorted(lam_hi):
    d = abs(mp.mpf(lam_hi[lbl]) - mp.mpf(lam_old[lbl]))
    print(lbl, "|new-old| =", mp.nstr(d, 2))
pickle.dump(lam_hi, open("lam_hi.pkl","wb"))
print("saved lam_hi.pkl (460 digits)")
