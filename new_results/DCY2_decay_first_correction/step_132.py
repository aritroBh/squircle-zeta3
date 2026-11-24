import os
base = "<project>/"
p = base+"research/08_theory/DCY2R_referee_varpi1/result.md"
print(os.path.exists(p))
print(open(p).read()[:9000])
