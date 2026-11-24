import os
base = '<project>/research/08_theory'
print(sorted(os.listdir(base)))
wts = f'{base}/WTS_weight_k_structure'
if os.path.exists(wts):
    for root, dirs, files in os.walk(wts):
        level = root.replace(wts,'').count(os.sep)
        if level < 2:
            print(root.replace(base,''), sorted(files)[:20])
