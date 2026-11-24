import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, workdir); sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import tow2r_lib as RL
import buildintegrand as BI
import genintegrand as GI

S = GI.S_poly(4, 4)
St_mine, D1 = RL.Stilde(S)
St_ref, D2 = BI.Stilde(S)
print('D:', D1, D2)
print('equal:', St_mine == St_ref)
print('mine:', [str(c) for c in St_mine])
print('ref :', [str(c) for c in St_ref])
