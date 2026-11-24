import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
os.chdir(os.path.join(workdir, 'tow2_copy'))
from fractions import Fraction as Fr
import towlib
towlib.WTS = '<project>/research/08_theory/WTS_weight_k_structure/theory'
import genintegrand as GI
from fastred import QSI, ZERO, ONE, IQ

exec(open('step_54.py').read().split("bad = []; tested = 0")[0].split("import towlib")[1].split("def rho_first_principles")[0])  # reuse helper defs? too fragile
