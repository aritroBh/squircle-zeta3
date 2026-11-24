import os
base = '<project>/research/08_theory'
src = open(os.path.join(base,'COR13C_third_corner_certificate','verify_cor13c.py')).read()
i = src.find('from fractions import Fraction as Fq')
print(src[i:i+5200])
