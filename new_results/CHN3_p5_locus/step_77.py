import sys; sys.path.insert(0,'.')
import numpy as np
import modengine as ME, point_engine as PE
import kimilib
p=2147483647
a=(1,1,1,1,1)
a_arrays=tuple(np.array([a[i]%p],dtype=np.int64) for i in range(5))
# word-level Areg W=5
vw_scalar = PE.volreg_words_point(a,5)[5]
vw_mod = ME.volreg_mod(a_arrays,5,p)[5]
print("scalar Areg[5]:", {w:int(v) for w,v in vw_scalar.items()})
print("mod    Areg[5]:", {w:int(v[0]) for w,v in vw_mod.items()}, " (mod p; p-1 means -1)")
# m0 W=5
m0s = kimilib.m0_coords(a,5)[5]
m0m = ME.m0_coords_mod(a_arrays,5,p)[5]
print("m0 scalar:", [str(x) for x in m0s])
print("m0 mod   :", [int(x[0]) for x in m0m])
