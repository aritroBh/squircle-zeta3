import os, sys, time
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, "<project>/research/08_theory/M4c_w5_closed_form")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
t0 = time.time()
import w5_cobracket as CB
print("import w5_cobracket OK in %.1fs" % (time.time()-t0))
from w4_symbol import QS2 as Q
ZERO, ONE = Q(0), Q(1)
P_LET, Q_LET = Q(2,1), Q(-1,-1)

t0 = time.time()
d = CB.delta_word(ZERO, (P_LET, Q_LET, ZERO, ONE, P_LET), ONE)
print("delta_word(p,q,0,1,p) in %.1fs:" % (time.time()-t0))
print(d)
t0 = time.time()
d2 = CB.delta_word(ZERO, (ZERO, P_LET, Q_LET, ONE, P_LET), ONE)
print("delta_word(0,p,q,1,p) in %.1fs:" % (time.time()-t0))
print(d2)
