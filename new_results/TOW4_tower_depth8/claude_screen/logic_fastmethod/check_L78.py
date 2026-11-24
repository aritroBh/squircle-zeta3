import sys, pickle, time
sys.path.insert(0, 'lib')
import sympy as sp
import tow2lib as T2, tow4lib as T4
d, m = T2.d, T2.m
Lraw = pickle.load(open('lib/L_simplified.pkl','rb'))
for r in range(1, 9):
    t0=time.time()
    cache = sp.sympify(Lraw[r], locals={'d':d,'m':m,'I':sp.I})
    diff = sp.cancel(sp.together(cache - T2.L(r)))
    print(f'L cache r={r} == tow2lib.L(r) symbolic: {diff==0}  ({round(time.time()-t0,1)}s)', flush=True)
print('check_tau8:', T4.check_tau8())
# tau table n=0..8 against DL definition numerically (50 digits) and conj symmetry
import mpmath as mp
mp.mp.dps = 50
for n in range(0, 9):
    t1 = complex(0)
    v1 = sp.N(T2.tau1(n), 50); v2 = sp.N(T2.tau2(n), 50)
    w1 = sp.N(T4.tau1f(n), 50); w2 = sp.N(T4.tau2f(n), 50)
    ok = abs(v1-w1) < 1e-40 and abs(v2-w2) < 1e-40
    print('tau', n, 'table==definition (50 digits):', ok)
