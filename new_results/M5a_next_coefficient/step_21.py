import os
base = '<project>'
print(sorted(os.listdir(f'{base}/paper/squircle_note/outputs')))
work = f'{base}/research/08_theory/M5a_next_coefficient'
# copy the tails machinery into work dir
import shutil
for f in ['m3_tails_exact.py', 'polylog_exact.py', 'w2_integrals.py', 'kummer.py', 'dilog_exact.py']:
    shutil.copy(f'{base}/paper/squircle_note/scripts/{f}', f'{work}/{f}')
print("copied tails machinery")