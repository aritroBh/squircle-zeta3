import os
base = '<project>/research/08_theory'
src = open(os.path.join(base,'COR7_classical_corners','result.md')).read()
print(src[:7000])
