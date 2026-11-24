import json
from fractions import Fraction as F
d=json.load(open('../M4e_w5_lchi/lchi_final_1200.json'))
c={nm:[F(a),F(q)] for nm,a,q in zip(d['names'],d['ca'],d['cb'])}
B2,Bp,Bm='Li5((3-2*sqrt2))','Li5((-1+1*sqrt2))','Li5((1-1*sqrt2))'
a,q=c.pop(B2)
for k in (Bp,Bm): c[k][0]+=16*a; c[k][1]+=16*q
a,q=c.pop('Li5(-1)')
c['z5'][0]+=F(-15,16)*a; c['z5'][1]+=F(-15,16)*q
c={nm:v for nm,v in c.items() if v[0] or v[1]}
c['Lchi8_5']=[F(d['aL']),F(d['bL'])]
saved=json.load(open('../M4e_w5_lchi/planner/R5_nice_reduced.json'))['terms']
saved={k:[F(v[0]),F(v[1])] for k,v in saved.items()}
print('keys equal:', set(c)==set(saved))
for k in sorted(set(c)|set(saved)):
    if c.get(k)!=saved.get(k):
        print('DIFF',k,'fresh:',c.get(k),'saved:',saved.get(k))
