base = '<project>/research/08_theory'
work = f'{base}/TOW3_tower_depth6'
import os
print(os.path.exists(f'{work}/progress.log'))
for f in ['step_02.py','step_03.py','step_04.py','step_05.py','step_06.py','step_07.py']:
    print('='*70); print(f, os.path.getsize(f'{work}/{f}'))
    print(open(f'{work}/{f}').read()[:2000])
