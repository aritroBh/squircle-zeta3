import os
base = '<project>/research/08_theory'
for f in ['WTS_weight_k_structure/theory/analyze_unfolded.out']:
    print(open(f'{base}/{f}').read()[:3000])
print('='*40)
# does the unfolded tower data exist?
print([f for f in os.listdir(f'{base}/WTS_weight_k_structure/theory') if 'unfold' in f])
