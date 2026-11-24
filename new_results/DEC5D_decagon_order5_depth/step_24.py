
import json
base = "<project>/research/08_theory/M4c_w5_closed_form/"
print(open(base+"symbol.json").read()[:1500])
print("="*60)
print(open(base+"cobracket.json").read()[:2000])
print("="*60)
import os
print(os.listdir(base+"planner"))
