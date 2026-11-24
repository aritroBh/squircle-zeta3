import sympy as sp, time, pickle, sys
t0=time.time()
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
import tow2lib as T2
d, m = T2.d, T2.m
i = sp.I
w = sp.Symbol('w')
def parity_parts(expr):
    ee = sp.expand(expr.subs({sp.sqrt(2): w}))
    P = sp.Poly(ee, w)
    a = sp.Integer(0); b = sp.Integer(0)
    for (e2,), c in P.terms():
        if e2 % 2 == 0: a += c * 2**(e2//2)
        else: b += c * 2**((e2-1)//2)
    return a, b
def reduce_Qi_const(expr):
    g = sp.together(expr)
    num, den = sp.fraction(g)
    na, nb = parity_parts(num)
    da, db = parity_parts(den)
    s2coef = sp.expand(nb*da - na*db)
    assert sp.expand(s2coef)==0, f"sqrt2 coeff nonzero"
    val = sp.cancel(sp.expand(na*da - 2*nb*db) / sp.expand(da**2 - 2*db**2))
    return sp.expand(val)

tau_tab = {}
for n in [8, 9, 10]:
    t1 = reduce_Qi_const(T2.tau1(n))
    t2 = reduce_Qi_const(T2.tau2(n))
    tau_tab[n] = (t1, t2)
    conj_ok = sp.expand(sp.conjugate(t1) - (-1)**n * t2) == 0
    print(f"tau1({n}) = {t1}   tau2({n}) = {t2}   conj-sym: {conj_ok}")
import mpmath as mp
mp.mp.dps = 40
def tau_num(which, n):
    s = n+1
    i_ = mp.mpc(0,1); a_ = mp.sqrt(2)+1; b_ = mp.sqrt(2)-1
    if which==1: return (-1)**s * ( (i_-b_)**(-(s+1)) + (i_+a_)**(-(s+1)) - i_**(-(s+1)) )
    return (-1)**s * ( (i_-a_)**(-(s+1)) + (i_+b_)**(-(s+1)) - i_**(-(s+1)) )
for n in [8,9,10]:
    for which in [1,2]:
        vn = tau_num(which, n)
        ve = tau_tab[n][which-1]
        diff = abs(mp.mpc(float(sp.re(ve)), float(sp.im(ve))) - vn)
        print(f"numeric tau{which}({n}): |diff| = {mp.nstr(diff, 5)}")
pickle.dump({n: (str(tau_tab[n][0]), str(tau_tab[n][1])) for n in tau_tab}, open(work+'/tau910.pkl','wb'))
print("total", round(time.time()-t0,1))
