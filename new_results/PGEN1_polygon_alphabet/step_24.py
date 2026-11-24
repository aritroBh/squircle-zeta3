import os
base = "<project>"
t = open(os.path.join(base,"research/08_theory/DEC1_decagon/result.md")).read()
print(t[:5500])
