base = '<project>'
m4b = base + '/research/08_theory/M4b_w4pp_exact'
for f in ['step03_numeric_checks.py','step04_stage1_modp.py']:
    print('#'*20, f, '#'*20)
    print(open(m4b+'/'+f).read())
