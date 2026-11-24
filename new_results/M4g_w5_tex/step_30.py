import re, json
from fractions import Fraction
# reproduce parse_atom on the suspicious line
line=open('supplement_W5.tex').read()
for l in line.splitlines():
    if '3061997' in l: print('LINE:', l)
# test split
body=r'-\bigl(\tfrac{3061997}{483840}+\tfrac{208211}{967680}\sqrt2\bigr)\,\lambda\,\pi^{2}\,\pi^{2}'
coeff,atom=body.split('\\,',1)
print('COEFF:',coeff)
print('ATOM :',atom)
print('PARTS:',atom.split('\\,'))
