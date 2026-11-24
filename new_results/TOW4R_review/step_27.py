
import os, re, pickle
import sympy as sp
work = "<project>/research/08_theory/TOW4R_review"
tow4 = "<project>/research/08_theory/TOW4_tower_depth8"
txt = open(tow4 + "/proof_depth8.md").read()
seg = txt[txt.index("a_7 = N_7"):txt.index("a_8 = N_8")]
print(seg)
