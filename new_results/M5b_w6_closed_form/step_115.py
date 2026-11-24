import time, os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
time.sleep(780)
print(open(m5b+'/w6_lll_1900.out').read()[-500:])
print('json1900:', os.path.exists(m5b+'/w6_lll_final_1900.json'))
