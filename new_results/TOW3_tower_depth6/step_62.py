base = '<project>/research/08_theory'
import os
wt5 = f'{base}/WT5_weight5_all_k'
print(sorted(os.listdir(wt5))[:30])
src = open(f'{wt5}/result.md').read()
print(src[:3500])
