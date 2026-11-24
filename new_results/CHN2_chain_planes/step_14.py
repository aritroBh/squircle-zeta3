
import pickle
from fractions import Fraction as F
chn1 = '<project>/research/08_theory/CHN1_chain_bodies'
forms = pickle.load(open(f'{chn1}/forms.pkl','rb'))
result = pickle.load(open('partB_planecoords.pkl','rb'))

def restrict4(poly4, plane):
    out={}
    for (e1,e2,e3,e4),c in poly4.items():
        k=(e1+e2, e3+e4) if plane=='Pi1' else (e1+e3, e2+e4)
        out[k]=out.get(k,F(0))+c
    return {k:c for k,c in out.items() if c}

print(f"{'W':>3} {'plane':>5} {'#basis':>6} {'maxdiff':>10}")
allok=True
for W in range(2,12):
    for plane in ['Pi1','Pi2']:
        d=len(result[plane][W])
        maxc=F(0); nb=0
        for b in range(d):
            f4=forms[(W,b)]
            fr=restrict4(f4, plane)
            mine=result[plane][W][b]
            diff=dict(fr)
            for e,c in mine.items():
                diff[e]=diff.get(e,F(0))-c
                if diff[e]==0: del diff[e]
            if diff:
                nb+=1
                for c in diff.values():
                    if abs(c)>maxc: maxc=c
        status = 'OK' if nb==0 else f'MISMATCH in {nb} basis elems, max {maxc}'
        if nb: allok=False
        print(f"{W:>3} {plane:>5} {d:>6}  {status}")
print("ALL EXACT MATCH W=2..11:", allok)
