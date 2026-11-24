import re
src = open('cor8r_sym.py').read()
i = src.find('def factor_Qc')
print(src[i:i+2200])
