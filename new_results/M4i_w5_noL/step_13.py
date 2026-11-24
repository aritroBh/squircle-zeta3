base = '<project>'
m4f = base + '/research/08_theory/M4f_w5_extras'
src = open(m4f + '/w5_search.py').read()
i = src.find('def evaluate')
print(src[i:i+2500])
print('='*20, 'w5_cobracket.py head')
c = open(m4f + '/w5_cobracket.py').read()
print(c[:2500])
