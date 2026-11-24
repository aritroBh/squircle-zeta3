import os
base = '<project>'
with open(os.path.join(base, 'paper/ieee/main.tex')) as f:
    lines = f.readlines()
# Find section IX "How the Proofs Work"
for i, l in enumerate(lines):
    if 'How the Proofs Work' in l or 'weight-two' in l.lower() or 'Rw2' in l:
        print(i+1, l.rstrip())
