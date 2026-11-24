import os, sys, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
for p in ["<project>/research/08_theory/M4c_w5_closed_form",
          "<project>/research/08_theory/A3_W4_proof",
          "<project>/paper/squircle_note/scripts"]:
    sys.path.insert(0, p)
import w5_cobracket as CB
from qsqrt2 import QS2 as Q
ZERO, ONE = Q(0), Q(1)
P_LET, Q_LET = Q(2,1), Q(-1,-1)

t0 = time.time()
d = CB.delta_word(ZERO, (P_LET, Q_LET, ZERO, ONE, P_LET), ONE)
print("delta_word(p,q,0,1,p) in %.1fs:" % (time.time()-t0))
print(d)
