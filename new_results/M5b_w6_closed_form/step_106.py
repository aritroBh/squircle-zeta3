import time, os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
time.sleep(900)
print(open(m5b+'/w6_lll_1600.out').read()[-700:])
print('json1600:', os.path.exists(m5b+'/w6_lll_final_1600.json'))
