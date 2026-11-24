base = '<project>/research/08_theory'
d = base + '/ZQ5d_zagier_d13_clean'
for f in ['step_01.py','step_02.py','step_03.py']:
    print('=====', f)
    print(open(d+'/'+f).read()[:2500])
