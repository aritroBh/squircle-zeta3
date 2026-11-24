import time, os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
time.sleep(240)
print(open(m5b+'/w6_lll_1300.out').read()[-400:])
print('json1300:', os.path.exists(m5b+'/w6_lll_final_1300.json'))
print('\n'.join(open(m5b+'/progress.log').read().splitlines()[-3:]))
