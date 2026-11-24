import sys; sys.path.insert(0,'.')
import numpy as np, importlib
import modengine as ME
importlib.reload(ME)
p=2147483647
pts=[(1,1,1,1,1)]
a_arrays=tuple(np.array([a[i]%p for a in pts],dtype=np.int64) for i in range(5))
res=ME.full_coords_mod(a_arrays,3,p)
print("W=3 z3 mod p:", int(res[3][0][0]), " expected 14")
# check volreg part alone
vw=ME.volreg_mod(a_arrays,3,p)
print("volreg W=3:", {w:int(v[0]) for w,v in vw[3].items()})
# expected Areg[3] at (1,1,1,1,1) = {(3,):-1} -> mod p = p-1
print("p-1 =", p-1)
# m0 part
m0=ME.m0_coords_mod(a_arrays,3,p)
print("m0 W=3:", int(m0[3][0][0]), " expected 16")
