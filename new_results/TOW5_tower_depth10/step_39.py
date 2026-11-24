import sympy as sp, time, pickle, sys
t0=time.time()
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
import tow2lib as T2
d, m = T2.d, T2.m
i = sp.I
aq = T2.aq; bq = T2.bq

def reduce_Qi_const(expr):
    """reduce a constant in Q(sqrt2,i) to a+b*sqrt2 form, assert b==0, return a in Q(i)"""
    e = sp.expand(expr)
    w = sp.Symbol('w')
    ee = sp.expand(e.subs({sp.sqrt(2): w}))
    P = sp.Poly(ee, w)
    a = sp.Integer(0); b = sp.Integer(0)
    for (ee2,), c in P.terms():
        if ee2 % 2 == 0: a += c * 2**(ee2//2)
        else: b += c * 2**((ee2-1)//2)
    assert sp.expand(b) == 0, f"sqrt2 part nonzero: {b}"
    return sp.expand(a)

tau_tab = {}
for n in [8, 9, 10]:
    t1 = reduce_Qi_const(T2.tau1(n))
    t2 = reduce_Qi_const(T2.tau2(n))
    tau_tab[n] = (t1, t2)
    # verify conjugation symmetry conj(tau1_n) = (-1)^n tau2_n
    conj_ok = sp.expand(sp.conjugate(t1) - (-1)**n * t2) == 0
    print(f"tau1({n}) = {t1}   tau2({n}) = {t2}   conj-sym: {conj_ok}")
# also 40-digit numeric cross-check vs tow2lib defs
import mpmath as mp
mp.mp.dps = 40
def tau1_num(n):
    s = n+1
    i_ = mp.mpc(0,1); a_ = mp.sqrt(2)+1; b_ = mp.sqrt(2)-1
    return (-1)**s * ( (i_-b_)**(-(s+1)) + (i_+a_)**(-(s+1)) - i_**(-(s+1)) )
for n in [9,10]:
    vn = tau1_num(n)
    ve = tau_tab[n][0]
    diff = abs(mp.mpc(float(sp.re(ve)), float(sp.im(ve))) - vn)
    print(f"numeric tau1({n}): |diff| = {mp.nstr(diff, 5)}")
pickle.dump({n: (str(tau_tab[n][0]), str(tau_tab[n][1])) for n in tau_tab}, open(work+'/tau910.pkl','wb'))
print("total", round(time.time()-t0,1))
