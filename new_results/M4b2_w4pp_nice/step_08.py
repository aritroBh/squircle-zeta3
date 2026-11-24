base = '<project>'
m4b = base + '/research/08_theory/M4b_w4pp_exact'
for f in ['step07_make_certificate.py','step08_numeric_sanity.py','verify_W4pp.py']:
    print('#'*20, f, '#'*20)
    print(open(m4b+'/'+f).read())
