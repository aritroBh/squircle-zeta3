base = '<project>'
m4f = base + '/research/08_theory/M4f_w5_extras'
src = open(m4f + '/subset_exact.py').read()
# print from 'def acol' onward
i = src.find('    def acol')
print(src[i:])
