import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 80
import numpy as np
av = 1+mp.sqrt(2); bv = mp.sqrt(2)-1
def tfun(v): return (v-bv)*(v+av)/(4*v)
f = lambda v: mp.log(tfun(v))/(v+av)
# slices toward v=b
tot = mp.mpf(0); x = av-bv
while x > mp.mpf(10)**(-68):
    tot += mp.quad(f, [bv+x/2, bv+x]); x /= 2
print("N(-a) slices:", mp.nstr(tot,15))
xs = np.linspace(float(bv)+1e-9, float(av)-1e-9, 2000000)
tf = (xs-float(bv))*(xs+float(av))/(4*xs)
est = np.trapezoid(np.log(tf)/(xs+float(av)), xs)
print("numpy estimate:", est)
# try plain quad on [b, a] with midpoint split
v2 = mp.quad(f, [bv, (av+bv)/2, av])
print("plain quad [b,mid,a]:", mp.nstr(v2,15))
