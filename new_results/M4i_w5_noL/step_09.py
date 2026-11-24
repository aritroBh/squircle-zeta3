base = '<project>'
m4f = base + '/research/08_theory/M4f_w5_extras'
src = open(m4f + '/optimize6.py').read()
i = src.find("log(f'seed-descended")
print(src[i:i+2200])
print('='*30, 'export_forms.py', '='*30)
print(open(m4f + '/export_forms.py').read())
