import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
for n in range(0, 8):
    t1 = sp.radsimp(sp.simplify(T.tau1(n)))
    t2 = sp.radsimp(sp.simplify(T.tau2(n)))
    print(f'tau1_{n} = {t1}    tau2_{n} = {t2}   in Q(i): {not t1.has(sp.sqrt(2)) and not t2.has(sp.sqrt(2))}')
