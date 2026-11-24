base = '<project>'
m4f = base + '/research/08_theory/M4f_w5_extras'
for f in ['m4f_common.py','overcomplete.py','drop_L.py']:
    print('='*30, f, '='*30)
    print(open(m4f + '/' + f).read())
