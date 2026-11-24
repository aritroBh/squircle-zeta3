import os
d = "<project>/research/08_theory/E5_edge_law/"
files = os.listdir(d)
print(files)
# find a result/readme
for f in files:
    if f.endswith(".md"):
        print("="*80); print(f)
        print(open(d+f).read()[:3000])
