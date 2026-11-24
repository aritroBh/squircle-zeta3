import os
base = "<project>"
t = open(os.path.join(base,"research/08_theory/DOD5_dodecagon_termwise/result.md")).read()
print(t[:7000])
