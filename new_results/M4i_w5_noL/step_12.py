base = '<project>'
m4f = base + '/research/08_theory/M4f_w5_extras'
src = open(m4f + '/w5_search.py').read()
print(src[:3000])
print('...TAIL...')
print(src[-1200:])
