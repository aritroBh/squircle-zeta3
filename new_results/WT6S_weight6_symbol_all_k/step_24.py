
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
# Check the head of local w4_symbol.py
src = open(os.path.join(WD, 'w4_symbol.py')).read()
print(src[:800])
print("...")
# Check qsqrt2 exports
q = open(os.path.join(WD, 'qsqrt2.py')).read()
import re
print("qsqrt2 defines:", re.findall(r'^(?:def |class )(\w+)', q, re.M)[:20])
print("A, B, R:", re.findall(r'^(A|B|R|SQ2)\s*=', q, re.M))
