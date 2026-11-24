
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
from fractions import Fraction as Fr
import reflib as R
import hexc2_common as HC

mine = R.my_li_atom(2, Fr(1,3))
theirs = HC.li_atom(2, Fr(1,3))
print("Li(2,1/3) equal:", mine == theirs, len(mine), len(theirs))
frame = HC.cross_table()[Fr(1,3)][0]
print("their rep:", frame)
print("my rep[0]:", R.cross_reps(Fr(1,3))[0])
frame_mine = tuple(None if p is HC.INF else p for p in frame)
mine2 = R.my_li_atom(2, Fr(1,3), frame=frame_mine)
print("same-frame equal:", mine2 == theirs)
