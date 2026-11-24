import os
base = '<project>/research/08_theory'
c7 = os.path.join(base,'COR7_classical_corners')
print(sorted(os.listdir(c7)))
print()
print(open(os.path.join(base,'COR13C_third_corner_certificate','cyclo.py')).read())
