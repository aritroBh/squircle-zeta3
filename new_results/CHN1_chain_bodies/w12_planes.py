
"""Weight-12 plane test: for points on Pi1 (s,s,1,1) and Pi2 (s,1,s,1), s=1..15,
compute c12 value coords in the W=12 basis, and the two quotient obstructions:
  obs1 = z39 + (14/9)*z57   (coefficient of z39 in MZV_12/span(products), using z57==(14/9)z39 mod products)
  obs2 = w*
Also membership rank test. Saves to w12_planes.json."""
import json, time
from fractions import Fraction as F
import mzv_exact as M, engine as E
def table(W):
    if W not in E._TABLES:
        E._TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return E._TABLES[W]
E.table=table
out=[]
for plane in ['Pi1','Pi2']:
    for s in range(1,16):
        a=(s,s,1,1) if plane=='Pi1' else (s,1,s,1)
        t=time.time()
        r=E.reduced_coeffs(a,12)[12]
        obs1=r[9]+F(14,9)*r[10]
        out.append({'plane':plane,'s':s,'z39':str(r[9]),'z57':str(r[10]),'w*':str(r[11]),
                    'z2^2z35':str(r[7]),'z2z37':str(r[8]),'obs1':str(obs1),'sec':round(time.time()-t,1)})
        json.dump(out, open('w12_planes.json','w'))
print('done', len(out))
