import sympy as sp, time, pickle, sys
t0=time.time()
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
import tow2lib as T2
d, m = T2.d, T2.m
i = sp.I
raw = pickle.load(open(work+'/R5_raw.pkl','rb'))
num5 = sp.sympify(raw['num5'], locals={'d':d,'m':m,'I':i})
den5 = sp.sympify(raw['den5'], locals={'d':d,'m':m,'I':i})

# sl_5 (proved = (log S)^(5)(t*))
R = T2.R
sl5 = 24*R[1]**5 - 60*R[1]**3*R[2] + 20*R[1]**2*R[3] + 30*R[1]*R[2]**2 - 5*R[1]*R[4] - 10*R[2]*R[3] + num5/den5
sl5c = sp.cancel(sl5)
print("sl5 cancelled, den factored:", sp.factor(sp.fraction(sl5c)[1])[:120] if len(str(sp.factor(sp.fraction(sl5c)[1])))>120 else sp.factor(sp.fraction(sl5c)[1]), f"({round(time.time()-t0,1)}s)")

# extend sl in tow2lib namespace and compute L(9), L(10)
T2.sl[5] = sl5c
t9=time.time()
L9 = T2.Lr(9)
print("L(9) built", round(time.time()-t9,1),"s")
t10=time.time()
L10 = T2.Lr(10)
print("L(10) built", round(time.time()-t10,1),"s")
# denominators
for nm, Lr in [('L9',L9),('L10',L10)]:
    nn, dd = sp.fraction(sp.together(Lr))
    print(nm, "den factored:", sp.factor(dd))
print("total", round(time.time()-t0,1))
