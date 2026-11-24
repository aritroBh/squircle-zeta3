import os
d = "<project>/research/08_theory/K1_nonlocal_law/"
for i in [1,2,3]:
    fn = d + f"step_{i:02d}.py"
    with open(fn) as f:
        print("="*30, fn)
        print(f.read()[:2000])
