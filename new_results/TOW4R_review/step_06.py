
import os
base = "<project>"
tow4 = base + "/research/08_theory/TOW4_tower_depth8"
txt = open(tow4 + "/result.md").read()
print(txt[3900:5400])
print("=====PROOF_DEPTH8=====")
print(open(tow4 + "/proof_depth8.md").read())
