
import pickle
from fractions import Fraction as F
result = pickle.load(open('partB_planecoords.pkl','rb'))
BASIS12 = ['z12','z3z9','z2z3z7','z2z5^2','z2^2z3z5','z2^3z3^2','z3^4','z2^2z(3,5)','z2z(3,7)','z(3,9)','z(5,7)','z(2,1,1,8)']
def fmt(poly):
    if not poly: return "0"
    terms=[]
    for (e1,e2),c in sorted(poly.items()):
        terms.append(f"{'+' if c>0 else ''}{c} x^{e1}y^{e2}")
    return ' '.join(terms)
for plane in ['Pi1','Pi2']:
    print(f"--- {plane} weight-12 coordinates ---")
    for b in [7,8,9,10,11]:
        print(f"  {BASIS12[b]:>12}: {fmt(result[plane][12][b])[:120]}")
    # check ratio r9 = -14/9 r10
    p9=result[plane][12][9]; p10=result[plane][12][10]
    # ratio test: p9 + 14/9 p10 == 0?
    combo={e:c+F(14,9)*p10.get(e,F(0)) for e,c in p9.items()}
    combo={e:c for e,c in combo.items() if c}
    for e,c in p10.items():
        if e not in p9: combo[e]=combo.get(e,F(0))+F(14,9)*c
    combo={e:c for e,c in combo.items() if c}
    print("  r9 + (14/9) r10 == 0 :", not combo)
