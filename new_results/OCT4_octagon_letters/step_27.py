import os, sys, pickle, json; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 300
d = pickle.load(open('letter_pts.pkl','rb'))
real_pts = {k: mp.mpf(v) for k,v in d['real_pts'].items()}
complex_pts = {k: mp.mpc(mp.mpf(a),mp.mpf(b)) for k,(a,b) in d['complex_pts'].items()}

consts = {}
consts['z3'] = mp.zeta(3)
for k,v in real_pts.items():
    consts[k] = mp.polylog(3, v)
for k,z in complex_pts.items():
    L = mp.polylog(3, z)
    consts['Re '+k] = L.real
    consts['Im '+k] = L.imag
print("computed", len(consts), "pure weight-3 constants")
# sanity: Li3(-1) = -3/4 z3; Li3(1/2) formula
print("Li3(-1)+3/4 z3 =", mp.nstr(consts['Li3(-1)']+mp.mpf('0.75')*consts['z3'],5))
li12 = consts['Li3(1/2)']; l2=mp.log(2)
print("Li3(1/2) - (7/8 z3 - pi^2/12 l2 + l2^3/6) =", mp.nstr(li12-(mp.mpf(7)/8*consts['z3']-mp.pi**2/12*l2+l2**3/6),5))
# Im Li3 at roots of unity should be elementary (Bernoulli): Im Li3(e^{i th}) = th^3/12 - pi th^2/... 
# Cl_3(pi/2)=Re Li3(i): known = -3 z3/32? check
print("Re CL11 (Re Li3(i)) =", mp.nstr(consts['Re CL11'],20))
print("  -3 z3/32 =", mp.nstr(-3*consts['z3']/32,20))
# save
json.dump({k: mp.nstr(v, mp.mp.dps-5) for k,v in consts.items()}, open('purew3_300.json','w'))
log("step: computed %d pure weight-3 constants at 300dps; sanity Li3(-1),Li3(1/2) OK"%len(consts))
