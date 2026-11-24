import os
tags = {11:'118900|worst_id', 12:'NS, NY = 64, 64', 13:'mean-value|dtruth', 14:'saddle containment', 15:'256x256|256', 16:'certify_H3.py output|d1', 17:'estimated full', 18:'supC = 0', 19:'verbatim', 20:'TIGHT box vs CODED'}
for i in range(1,22):
    fn = f'step_{i:02d}.py'
    src = open(fn).read()
    lab = []
    for k,v in tags.items():
        for pat in v.split('|'):
            if pat in src: lab.append(k); break
    first = [l for l in src.split('\n') if l.strip() and not l.startswith('import') and not l.startswith('sys')][:1]
    print(f"{fn}: match={lab}  first={first}")
