import sys, os
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, tow2); os.chdir(tow2)
import sympy as sp
sys.setrecursionlimit(200)
import tow2lib as T
try:
    e = T.e_desc(4,1)
    print('e_desc(4,1) OK:', sp.simplify(e) != None)
except RecursionError as ex:
    print('RecursionError — e_desc(4,1) does NOT terminate with current tow2lib')
# try e_desc(3,0) and (3,1)
for ds in [(2,0),(3,0),(3,1),(4,0),(4,2)]:
    try:
        sys.setrecursionlimit(400)
        e = T.e_desc(*ds)
        print(ds, 'ok')
    except RecursionError:
        print(ds, 'RecursionError')
