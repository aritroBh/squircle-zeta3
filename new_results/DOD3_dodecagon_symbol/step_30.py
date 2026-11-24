import sys, os, pickle, time
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 620
import sympy as sp
from dod3_R import *

base = os.path.dirname(os.getcwd())
d1 = os.path.join(base,'DOD1_dodecagon'); d2 = os.path.join(base,'DOD2_dodecagon_finish')
atoms620 = pickle.load(open(os.path.join(d1,'atoms620.pkl'),'rb'))
g = pickle.load(open(os.path.join(d2,'gdswaps620.pkl'),'rb'))

vals = {l2: mp.log(2), l3: mp.log(3), lu: mp.log(2-mp.sqrt(3))}
def ev(e):
    e = e.subs(vals)
    e = e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,700))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,700)))))
    return mp.mpf(sp.N(e, 620))

t0 = time.time()
Rs = {}; Rv = {}; Rq = {}
worst = mp.mpf(0)
for (i,j,k) in TRIPLES:
    R = build_R(i,j,k)
    Rs[(i,j,k)] = R
    rv = ev(R)
    rq = atoms620['T3a'][(i,j,k)] - (g['SwG'][(i,j,k)]-g['SwD'][(i,j,k)])/atoms620['X'][i]
    Rv[(i,j,k)] = rv; Rq[(i,j,k)] = rq
    d = abs(rv-rq); worst = max(worst,d)
    print("R(%d;%d,%d): |sym-quad| = %s   (t=%.0fs)"%(i,j,k, mp.nstr(d,3), time.time()-t0))
print("\nWORST |R_sym - R_quad| over 12 triples:", mp.nstr(worst,3), " at 620 dps")
pickle.dump({'Rs':{k:str(v) for k,v in Rs.items()}}, open('R_all12_sym_str.pkl','wb'))
# store evaluated too
pickle.dump({'Rv':Rv,'Rq':Rq}, open('R_all12_vals620.pkl','wb'))
print("saved")
