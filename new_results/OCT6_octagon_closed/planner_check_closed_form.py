# Planner check: evaluate OCT6's 35-term closed form (table_Stot.pkl) with independently computed Li2/Li3 values
# and compare 48u*S with the independent 200-digit P'''(1) of D15b (key "4", entry 3).
import pickle, json, mpmath as mp, sympy as sp
mp.mp.dps=210
s2=mp.sqrt(2); u=mp.sqrt(2+s2)
pts={'x1':1-s2,'x2':1-s2/2,'k':s2-1,'d':2-s2,'r2h':s2/2,'mr2h':-s2/2}
C={'1':mp.mpf(1),'l2':mp.log(2),'lu':mp.log(1+s2),'z2':mp.pi**2/6,'z3':mp.zeta(3)}
for k,a in pts.items(): C['L2'+k]=mp.polylog(2,a); C['L3'+k]=mp.polylog(3,a)
def val(name):
    v=mp.mpf(1)
    for f in name.split('*'):
        f=f.strip(); p=1
        if '^' in f: f,p=f.split('^'); p=int(p)
        v*=C[f]**p
    return v
T=pickle.load(open('table_Stot.pkl','rb'))
S=mp.mpf(0)
for w,name,c in T:
    cf=mp.mpf(sp.N(sp.sympify(c),220)) if c!=0 else 0
    S+=cf*val(name)
P3=48*u*S
ref=mp.mpf(json.load(open('../D15b_polygon_law/vals_m2_8_200.json'))['4'][3])
print('terms',len(T),' 48uS =',mp.nstr(P3,30),' ref',mp.nstr(ref,30),' diff',mp.nstr(P3-ref,3))
