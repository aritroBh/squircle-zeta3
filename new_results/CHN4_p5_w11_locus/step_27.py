import pickle
from flint import fmpq_mat, fmpq
d=pickle.load(open('interp_points.pkl','rb')); sel=d['sel']; orbits=d['orbits']; mons7=d['mons7']
orbit_of={}
for i,o in enumerate(orbits):
    for e in o: orbit_of[e]=i
NORB=len(orbits)
def orbit_vals_int(a):
    vals=[0]*NORB
    for e in mons7:
        pw=1
        for x,ei in zip(a,e):
            if ei: pw*= x**ei
        vals[orbit_of[e]]+=pw
    return vals
data=[fmpq(v) for a in sel for v in orbit_vals_int(a)]
Mm=fmpq_mat(NORB,NORB,data)
print("exact rank over Q of SELECTED set:", Mm.rank(), "(need", NORB, ")")
