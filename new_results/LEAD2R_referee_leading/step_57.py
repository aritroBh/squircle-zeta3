import sys; sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 80
from math import factorial, comb
import pickle
from fractions import Fraction as Fr
from lead2r_lib import e_coeffs_logE, zigzag_ints, eps_m
E = e_coeffs_logE(14)
A = zigzag_ints(12)
av = 1+mp.sqrt(2); bv = mp.sqrt(2)-1
POLEV = {'-a':-av, '0':mp.mpf(0), '-b':-bv, 'a':av}
def tfun(v): return (v-bv)*(v+av)/(4*v)
def omfun(v): return (av-v)*(v+bv)/(4*v)
def Nint(i, j, r):
    # int_b^a dv/(v-r) log^i t log^j(1-t); endpoint logs at v=a (1-t) and v=b (t); r outside (b,a) except none
    f = lambda v: mp.log(tfun(v))**i*mp.log(omfun(v))**j/(v-r)
    tot = mp.mpf(0)
    # slices toward v=b (from above) and v=a (from below)
    x = av-bv
    while x > mp.mpf(10)**(-68):
        tot += mp.quad(f, [bv+x/2, bv+x]); x /= 2
    y = av-bv
    while y > mp.mpf(10)**(-68):
        tot += mp.quad(f, [av-y, av-y/2]); y /= 2
    return tot
res = pickle.load(open('lead2r_Dtable.pkl','rb'))
Nc = {}
def N(i,j,pole):
    key=(i,j,pole)
    if key not in Nc: Nc[key] = Nint(i,j,POLEV[pole])
    return Nc[key]
rows = []
worst = mp.mpf(0)
for w in range(2, 12):
    d = w-1
    tab, spare = res[d]
    tot = mp.mpf(0)
    for (key,pole),(D,poly,s2,pred) in tab.items():
        if D:
            tot += mp.mpf(int(D.numerator))/int(D.denominator)*N(key[0],key[1],pole)
    cw = tot/factorial(d-1)
    law = eps_m(w-2)*A[w-1]/(2*factorial(w-2))*E[w]
    r = abs(cw-law); worst = max(worst,r)
    rows.append((w, d, spare, r))
    print(f"w={w:2d} (d={d:2d}, spare={spare}): |c(w)-law| = {mp.nstr(r,3)}")
print("worst:", mp.nstr(worst,3))
pickle.dump([(w, mp.nstr(r, 5)) for w,d,s,r in rows], open('lead2r_cw_resid.pkl','wb'))
with open('progress.log','a') as f:
    f.write("step_23: END-TO-END from raw kernels (own D-extraction incl d=10 spare=1, own N-quadrature dps80): c(w)==law for w=2..11, worst "+mp.nstr(worst,3)+"\n")
