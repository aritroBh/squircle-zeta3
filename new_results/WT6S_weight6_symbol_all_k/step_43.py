
import os, re
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
q = open(os.path.join(WD, 'qsqrt2.py')).read()
# print the arithmetic methods
for m in re.findall(r'def (__mul__|__rmul__|__add__|__radd__|__sub__|__neg__|__eq__|__init__|__hash__)\(self.*?(?=\n    def |\nclass |\Z)', q, re.S):
    print(m[0] if isinstance(m, tuple) else m)
    print('---')
# Just print first 60 lines of QS2 class
i = q.find('class QS2')
print(q[i:i+2500])
