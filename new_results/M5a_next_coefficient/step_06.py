base = '<project>'
p = f'{base}/research/08_theory/M5a_next_coefficient'
for f in ['step_01.py','step_02.py','step_03.py','step_04.py','step_05.py']:
    print('='*40, f)
    print(open(f'{p}/{f}').read())
