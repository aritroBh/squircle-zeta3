"""Side probe: weight-2 piece of M_1 (top weight of M_1) through the same pipeline -- is k = d = 1 an exception at weight 2?"""
import pickle
from w4_common import *
import buildintegrand
ipf = buildintegrand.build_integrand_pf(1)
pickle.dump(ipf, open(os.path.join(HERE, "binteg_1.pkl"), "wb"))
red = fastred.reduce_k_fast(1, integrand_pf=ipf)
pickle.dump(red, open(os.path.join(HERE, "red_1.pkl"), "wb"))
set_dps(90); print("mp.mp.dps =", mp.mp.dps)
from w4_known import W2_exact
for key in sorted(red): print(" key", key, "poles", sorted(str(r) for r in red[key][0]))
pcs = pieces_mp(1, 1)
val, err = integrate(pcs)
print("W_2(1) =", mp.nstr(val, 30), " (pi^2+6 log^2 2)/24 =", mp.nstr(W2_exact(), 30))
print("W_2(1) - law value =", mp.nstr(val - W2_exact(), 8), " dps", mp.mp.dps)
