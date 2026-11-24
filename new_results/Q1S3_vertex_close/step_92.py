import pickle, mpmath as mp
from q1s_exact import femb, fconj, feq
mp.mp.dps=60
sol = pickle.load(open('q3_solve.pkl','rb'))
sols = pickle.load(open('sunit_li3.pkl','rb'))
def ReLi3(i):
    z = femb(sols[i]['z'],1)
    return mp.re(mp.polylog(3,z))
# the generator value is 2*Re Li3(z_i) for pair (i,j)
tot = mp.mpf(0)
print("V3' ReLi3 combination:")
for (i,j),c in sol['V3p']['li3']:
    v = 2*ReLi3(i)
    tot += mp.mpf(float(c))*v
    print("  pair (%3d,%3s) c=%6s  z=%s  2ReLi3=%s"%(i,j,c,
          mp.nstr(femb(sols[i]['z'],1),8), mp.nstr(v,12)))
import json
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
ref3=mp.mpf(json.load(open(b3+"V3p_1250.json"))['V3p'])
print("sum c*2ReLi3 =", mp.nstr(tot,20))
print("ref V3'      =", mp.nstr(ref3,20))
print("difference   =", mp.nstr(ref3-tot,12))
