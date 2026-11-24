import pickle
from fractions import Fraction as F
qf=pickle.load(open('Q11_fit.pkl','rb')); coefs=qf['coefs']; orbits=qf['orbits']; mons7=qf['mons7']
ev=pickle.load(open('f11_evals.pkl','rb'))
# Build explicit Q11 polynomial: exponent -> Fraction
Q11={}
for i,o in enumerate(orbits):
    c=coefs[i]
    if c==0: continue
    for e in o:
        Q11[e]=Q11.get(e,F(0))+c
print("Q11 explicit terms:", len(Q11))
def peval(poly,a):
    a=[F(x) for x in a]; tot=F(0)
    for e,c in poly.items():
        pw=F(1)
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        tot+=c*pw
    return tot
def Lval(a):
    a1,a2,a3,a4,a5=a
    return a2*a3*a4*(a1+a2+a3+a4+a5)

# VERIFY on 48 held-out points
verpts=[p for p,(lb,v) in ev.items() if lb=='verify']
print("verification points:", len(verpts))
bad=0; maxerr=None
for a in verpts:
    lhs=ev[a][1]           # f11(a) exact
    rhs=Lval(a)*peval(Q11,a)
    if lhs!=rhs:
        bad+=1; maxerr=(a,lhs,rhs)
print("MISMATCHES:", bad)
if bad: print("example:", maxerr)
else:
    print("ALL 48 VERIFICATION POINTS MATCH EXACTLY: f11 = a2*a3*a4*(Sum a)*Q11")
# also re-verify on the two planner V10 points (approx) later; first save Q11
pickle.dump(Q11, open('Q11_poly.pkl','wb'))
print("saved Q11_poly.pkl")
