import os
d = "<project>/research/08_theory/K1_nonlocal_law/"
for i in [4,5,6,7,8]:
    fn = d + f"step_{i:02d}.py"
    with open(fn) as f:
        print("="*30, fn)
        print(f.read()[:1800])
