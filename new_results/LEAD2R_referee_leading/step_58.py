import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 80
import pickle
res = pickle.load(open('lead2r_Dtable.pkl','rb'))
tab, spare = res[1]
print("d=1 extracted table:")
for k, v in tab.items(): print("  ", k, "Dtop=", v[0], "pred=", v[3])
av = 1+mp.sqrt(2); bv = mp.sqrt(2)-1
def tfun(v): return (v-bv)*(v+av)/(4*v)
def Nint_dbg(i, j, r, label):
    f = lambda v: mp.log(tfun(v))**i/(v-r)
    tot = mp.mpf(0); x = av-bv
    n=0
    while x > mp.mpf(10)**(-68):
        tot += mp.quad(f, [bv+x/2, bv+x]); x /= 2; n+=1
    print(f"  N({label}) = {mp.nstr(tot,15)}")
    return tot
# test: int_b^a log t dv/(v+a) numerically vs crude midpoint estimate
import numpy as np
xs = np.linspace(float(bv)+1e-6, float(av)-1e-6, 200000)
tf = (xs-float(bv))*(xs+float(av))/(4*xs)
est = np.trapz(np.log(tf)/(xs+float(av)), xs)
print("crude numpy estimate N(-a):", est)
N = Nint_dbg(1,0,-av,'-a')
