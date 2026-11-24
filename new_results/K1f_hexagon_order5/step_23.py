
import sympy as sp, pickle, time
WD = '<project>/research/08_theory/K1f_hexagon_order5'
A0,A1,B,kap=sp.symbols('A0 A1 B kap'); L2,q0,V02=sp.symbols('L2 q0 V02')
lA0,lA1,lB=sp.symbols('lA0 lA1 lB'); t,X,Y,W=sp.symbols('t X Y W')
local={'A0':A0,'A1':A1,'B':B,'kap':kap,'L2':L2,'q0':q0,'V02':V02,'lA0':lA0,'lA1':lA1,'lB':lB,'t':t,'X':X,'Y':Y,'W':W}
Ea=pickle.load(open(f'{WD}/E5_general_approachA.pkl','rb'))['E']
t0=time.time()
polys={}
keyset=(t,X,Y,W)
for k in range(6):
    ek=sp.sympify(Ea[k],locals=local)
    D={}
    for term in sp.Add.make_args(sp.expand(ek)):
        pd=term.as_powers_dict()
        key=(pd.get(t,0),pd.get(X,0),pd.get(Y,0),pd.get(W,0))
        rest=[s**e for s,e in pd.items() if s not in keyset]
        co=sp.Mul(*rest) if rest else sp.Integer(1)
        if key in D: D[key]=D[key]+co
        else: D[key]=co
    D={kk:sp.expand(v) for kk,v in D.items()}
    polys[k]=D
    print(f"E[{k}]: {len(D)} monomials in (t,X,Y,W), {time.time()-t0:.0f}s", flush=True)
pickle.dump(polys,open(f'{WD}/E5_polys.pkl','wb'))
print("saved E5_polys.pkl")
# log
with open(f'{WD}/progress.log','a') as f:
    f.write("step_23: rebuilt E5_polys.pkl (monomial decomposition in (t,X,Y,W)); sizes %s\n" % {k:len(v) for k,v in polys.items()})
