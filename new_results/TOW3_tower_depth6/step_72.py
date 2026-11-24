base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
src = open(f'{wts}/fastred.py').read()
import re
for mm in re.finditer(r'DL[12][^\n]*\n', src):
    print(mm.group(0))
# find pf_mul_DL and pf_eval_at
i = src.find('def pf_mul_DL')
print(src[i:i+1500])
