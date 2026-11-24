import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, workdir); sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import tow2r_lib as RL
import towlib
towlib.WTS = '<project>/research/08_theory/WTS_weight_k_structure/theory'
mine = RL.rho_first_principles(4, 0, 4)
for delta in range(5):
    print(delta, 'mine:', mine[delta])
    print(delta, 'ref :', towlib.rho((4,0), delta))
