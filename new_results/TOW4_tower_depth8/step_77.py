import sys
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4, tow2lib as T2
print("zigzag A[1..7]:", [T2.A[n] for n in range(1,8)])
print("e_top(8) =", T2.e_top(8))
print("F(8,4) =", T4.forced(8,4))
print("e_top(7) =", T2.e_top(7), " F(7,3) =", T4.forced(7,3))
# fix the proof text
work = '<project>/research/08_theory/TOW4_tower_depth8'
pt = open(work + '/proof_depth8.md').read()
bad_start = pt.find("  F(8,4) = 14 / Pi_4")
bad_end = pt.find("[The tops are proved == zigzag")
if bad_start != -1:
    pt = pt[:bad_start] + "  F(8,4) = 16 / Pi_4   (zigzag top, (-1)^4 A_5/Pi_4 with A_5 = 16)\n" + pt[pt.find("[The tops"):]]
open(work + '/proof_depth8.md','w').write(pt)
print("fixed F(8,4) line")
