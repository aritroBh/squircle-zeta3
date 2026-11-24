import time, os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
time.sleep(600)
print(open(m5b+'/w6_lll_1300b.out').read()[-600:])
print('json:', os.path.exists(m5b+'/w6_lll_final_1300.json'))
