import os, glob
base = '<project>'
th = os.path.join(base, 'research/08_theory')
print([d for d in os.listdir(th) if os.path.isdir(os.path.join(th,d))])
# look for weight-4 / WT4
for d in os.listdir(th):
    if 'WT4' in d or 'weight4' in d.lower() or 'weight_4' in d.lower():
        print("FOUND", d, os.listdir(os.path.join(th,d)))
