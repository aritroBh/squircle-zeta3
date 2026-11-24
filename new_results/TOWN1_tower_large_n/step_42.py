import os
base = '<project>'
p = os.path.join(base, 'research/08_theory/WT4_weight4_all_k')
print(os.listdir(p))
for fn in os.listdir(p):
    if fn.endswith('.md'):
        print("="*40, fn)
        with open(os.path.join(p, fn)) as f:
            print(f.read()[:6000])
