WORK='<project>/research/08_theory/M4f_w5_extras'
for f in ['step_01.py','step_02.py']:
    print('='*30, f)
    print(open(f'{WORK}/{f}').read()[:3000])
