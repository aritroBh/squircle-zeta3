import os
base = '<project>/research/08_theory'
src = open(os.path.join(base,'HEXC1_hexagon_certificate','cmzv6.py')).read()
print(src[:9500])
