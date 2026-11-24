import os
base = "<project>/"
p = base+"research/08_theory/DCY2R_referee_varpi1/result.md"
txt = open(p).read()
print(txt[9000:13000])
print("="*80)
import os
d = base+"research/08_theory/DCY2R_referee_varpi1/"
print(os.listdir(d))
