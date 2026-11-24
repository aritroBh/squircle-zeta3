
import os
base = "<project>"
tow4 = base + "/research/08_theory/TOW4_tower_depth8"
txt = open(tow4 + "/result.md").read()
print(txt[5200:9800])
