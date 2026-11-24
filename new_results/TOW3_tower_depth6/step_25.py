import sys, os, time
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
t0 = time.time()
try:
    e41 = T.e_desc(4, 1)
    print('e_desc(4,1) computed in', round(time.time()-t0,1), 's')
    print(sp.cancel(e41))
except RecursionError:
    print('RECURSION ERROR')
