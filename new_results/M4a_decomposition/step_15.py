import os
BASE = '<project>'
W4D = os.path.join(BASE,'research/08_theory/A3_W4_proof')
src = open(os.path.join(W4D,'w4_symbol.py')).read()
i = src.find('def w4_g_expression')
print(src[i:i+2200])
