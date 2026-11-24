import os, json
from fractions import Fraction as Ff
print('loadavg:', os.getloadavg())
J = json.load(open('out/d13/candidates.json'))
allc = J['cands']; names=J['names']
def xf(c): return abs(float(Ff(c['x'][0]) + Ff(c['x'][1])*3.6055512754639892931))
sortedc = sorted(allc, key=lambda c: (c['h'], xf(c)))
cap = json.load(open('out/d13/candidates_cap.json'))
# verify first 275 match cap-275 pool exactly
match = all(sortedc[i]['x']==cap['cands'][i]['x'] and sortedc[i]['vx']==cap['cands'][i]['vx'] for i in range(275))
print('first 275 of re-sorted full pool == cap275 pool:', match)
if not match:
    bad=[i for i in range(275) if sortedc[i]['x']!=cap['cands'][i]['x']]
    print('mismatch idx sample:', bad[:5])
    # maybe sort ties differ; check as sets
    s1={tuple(c['x']) for c in sortedc[:275]}; s2={tuple(c['x']) for c in cap['cands']}
    print('as sets equal:', s1==s2, 'sizes', len(s1), len(s2))
c350 = sortedc[:350]
json.dump(dict(D=13, B='cap', split=[3], inert=[2], names=names, S=[13,2,3], cap=350,
               cands=[dict(x=c['x'], vx=c['vx'], v1m=c['v1m'], h=c['h'], xf=xf(c)) for c in c350]),
          open('out/d13/candidates_cap350.json','w'))
from collections import Counter
print('cap350 saved; heights:', dict(sorted(Counter(c['h'] for c in c350).items())), 'neg x:', sum(1 for c in c350 if xf(c)!=0 and float(Ff(c['x'][0])+Ff(c['x'][1])*3.6055512754639) < 0))
