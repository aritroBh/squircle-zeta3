# s11: is the cached L_parsed.pkl (used for a_5, a_6) EXACTLY the proved construction tow2lib.L(r)? (symbolic, not random-point)
# and does Bell(L_parsed) = tow2lib.complete_bell?  tow2lib.L(r) is built from R_1..R_4 (proved by prove_S12/S34) + algebraic pieces.
import sys, time, pickle
sys.dont_write_bytecode = True
import sympy as sp
import tow2lib as T2, tow3lib as T3
d, m = T3.d, T3.m
t0 = time.time()
Lp = T3._load_L()
print('L_parsed keys:', sorted(Lp.keys()), type(Lp[1]))
ok = {}
for r in range(1, 7):
    ref = T2.L(r)
    diff = sp.together(sp.expand(ref - Lp[r]))
    num = sp.numer(diff)
    # radicals: simplify numerator coefficientwise
    num = sp.expand(sp.radsimp(sp.expand(num)))
    ok[r] = (num == 0) or (sp.simplify(num) == 0)
    print(f'  L({r}): L_parsed - tow2lib.L(r) == 0 : {ok[r]}   ({round(time.time()-t0,1)}s)', flush=True)
print('all r<=6 equal:', all(ok.values()))
for dl in (5, 6):
    ref = T2.arho(dl)      # complete Bell of tow2lib.L(1..dl)
    diff = sp.together(sp.expand(ref - T3.arho(dl)))
    num = sp.expand(sp.radsimp(sp.expand(sp.numer(diff))))
    print(f'  a_{dl}: Bell(tow2lib.L) - tow3lib.arho == 0 : {num == 0 or sp.simplify(num)==0}   ({round(time.time()-t0,1)}s)', flush=True)
