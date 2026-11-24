import os
print('CWD:', os.getcwd())
base = '<project>'
m4b = base + '/research/08_theory/M4b_w4pp_exact'
m4a = base + '/research/08_theory/M4a_decomposition'
print(open(m4b+'/m4b_common.py').read())
print('='*80)
print(open(m4a+'/w4pp_closed_form.json').read())
