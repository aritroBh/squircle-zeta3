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
R = T2.R
sl5 = sp.cancel(24*R[1]**5 - 60*R[1]**3*R[2] + 20*R[1]**2*R[3] + 30*R[1]*R[2]**2 - 5*R[1]*R[4] - 10*R[2]*R[3] + num5/den5)
T2.sl[5] = sl5
w = sp.Symbol('w')
def split_sqrt2(expr):
    g = sp.together(expr)
    num, den = sp.fraction(g)
    nume = sp.expand(num.subs({sp.sqrt(2): w}))
    Pn = sp.Poly(nume, w)
    a = sp.Integer(0); b = sp.Integer(0)
    for (e,), c in Pn.terms():
        if e % 2 == 0: a += c * 2**(e//2)
        else: b += c * 2**((e-1)//2)
    # denominator may also have sqrt2: rationalize
    dene = sp.expand(den.subs({sp.sqrt(2): w}))
    Pd = sp.Poly(dene, w)
    da = sp.Integer(0); db = sp.Integer(0)
    for (e,), c in Pd.terms():
        if e % 2 == 0: da += c * 2**(e//2)
        else: db += c * 2**((e-1)//2)
    return (a, b), (da, db)

def simplify_Qi(expr):
    """reduce expr in Q(sqrt2,i)(d,m) to Q(i)(d,m) if sqrt2-part vanishes"""
    (a,b),(da,db) = split_sqrt2(expr)
    assert sp.expand(b)==0, "numerator has sqrt2 part!"
    # value = a/(da+db*sqrt2) = a(da-db*sqrt2)/(da^2-2db^2); sqrt2 coeff must vanish
    den2 = sp.expand(da**2 - 2*db**2)
    num_a = sp.expand(a*da); num_b = sp.expand(-a*db)
    assert sp.expand(num_b)==0, "value still has sqrt2!"
    return sp.cancel(num_a/den2)

res = {}
for r in [9, 10]:
    t1=time.time()
    Lr = T2.Lr(r)
    s = simplify_Qi(Lr)
    nn, dd = sp.fraction(sp.together(s))
    print(f"L({r}): simplified to Q(i)(d,m); den = {sp.factor(dd)}  ({round(time.time()-t1,1)}s)")
    res[r] = s
pickle.dump({r: str(res[r]) for r in res}, open(work+'/L910_candidate.pkl','wb'))
print("saved; total", round(time.time()-t0,1))
