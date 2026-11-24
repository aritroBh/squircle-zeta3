base = '<project>/research/'
for f in ['cor6_nf.py', 'cor6_sunit.py', 'cor6_sym.py', 'cor6_ft.py', 'cor6_pipe.py', 'cor6_table.py']:
    print('='*30, f, '='*30)
    print(open(base+'08_theory/COR6_true_point_symbols/'+f).read())
