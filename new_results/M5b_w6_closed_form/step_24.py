import os, shutil
m5b='<project>/research/08_theory/M5b_w6_closed_form'
m4c='<project>/research/08_theory/M4c_w5_closed_form'
m4e='<project>/research/08_theory/M4e_w5_lchi'
m4f='<project>/research/08_theory/M4f_w5_extras'
m5a='<project>/research/08_theory/M5a_next_coefficient'
a3='<project>/research/08_theory/A3_W4_proof'
# copies (never edit originals)
for src,dst in [
    (f'{m5a}/omega5.json','omega5.json'),
    (f'{m5a}/W6_400.txt','W6_400.txt'),
    (f'{m5a}/quad5_100.json','quad5_100.json'),
    (f'{m5a}/n6_solution.json','n6_solution.json'),
    (f'{m4c}/w5_symbol.py','w5_symbol.py'),
    (f'{m4c}/w5_cobracket.py','w5_cobracket.py'),
    (f'{m4c}/w5_search.py','w5_search.py'),
    (f'{m4e}/nice_search.py','nice_search.py'),
    (f'{m4e}/Lchi8_5.json','Lchi8_5.json'),
    (f'{m4e}/cobracket.json','cobracket_M4.json'),
    (f'{a3}/w4_symbol.py','w4_symbol.py'),
    (f'{a3}/w4_right_coproduct.py','w4_right_coproduct.py'),
    (f'{m4f}/R5_min3.json','R5_min3.json'),
    (f'{m4f}/pool.json','pool.json'),
]:
    shutil.copy(src, os.path.join(m5b,dst))
    print('copied', dst)
# check w4_right_coproduct imports
head=open(f'{a3}/w4_right_coproduct.py').read()[:1500]
print(head)
