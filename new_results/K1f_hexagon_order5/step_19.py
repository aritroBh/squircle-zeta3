
import sympy as sp, pickle, time
WD = '<project>/research/08_theory/K1f_hexagon_order5'
A0,A1,B,kap=sp.symbols('A0 A1 B kap'); L2,q0,V02=sp.symbols('L2 q0 V02')
lA0,lA1,lB=sp.symbols('lA0 lA1 lB'); t,X,Y,W=sp.symbols('t X Y W')
local={'A0':A0,'A1':A1,'B':B,'kap':kap,'L2':L2,'q0':q0,'V02':V02,'lA0':lA0,'lA1':lA1,'lB':lB,'t':t,'X':X,'Y':Y,'W':W}
Ea=pickle.load(open(f'{WD}/E5_general_approachA.pkl','rb'))['E']
t0=time.time()
polys={}
for k in range(6):
    ek=sp.sympify(Ea[k],locals=local)
    P=sp.Poly(sp.expand(ek),t,X,Y,W)
    polys[k]={tuple(m):c for m,c in P.terms()}
    print(f"E[{k}]: {len(polys[k])} monomials, {time.time()-t0:.0f}s")
pickle.dump(polys,open(f'{WD}/E5_polys.pkl','wb'))
print("saved E5_polys.pkl")
