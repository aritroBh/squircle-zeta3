exec(open('common.py').read())
C4, Q7 = get_C4_Q7()
forms = load_forms()
P8d = forms[8][3]
revd = lambda p: {(e[4],e[3],e[2],e[1],e[0]):c for e,c in p.items()}
print("P8 == rev(P8):", revd(P8d)==P8d)
# check C4 invariance numerically
import random
for pt in [(1,2,3,4,5),(2,7,1,3,5),(1,1,2,1,3)]:
    rpt = pt[::-1]
    print(pt, "C4:", C4.subs(dict(zip(A,pt))), " rev:", C4.subs(dict(zip(A,rpt))),
          " Q7':", Q7.subs(dict(zip(A,pt))), " rev:", Q7.subs(dict(zip(A,rpt))))
