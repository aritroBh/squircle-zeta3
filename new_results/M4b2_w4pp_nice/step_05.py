base = '<project>'
m4b = base + '/research/08_theory/M4b_w4pp_exact'
for f in ['step_01.py','step02_build_atoms.py']:
    print('#'*20, f, '#'*20)
    print(open(m4b+'/'+f).read())
