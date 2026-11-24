import os, subprocess
base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
files = sorted(os.listdir(wts))
print([f for f in files if f.endswith('.md') or f.endswith('.txt')])
# grep for Lemma 3 in any md/py
out = subprocess.run(['grep','-rln','Lemma 3\|Lemma 2\|Lemma 4', wts], capture_output=True, text=True)
print(out.stdout[:2000])
