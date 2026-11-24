import os
for i in [4,5,6,7,8]:
    fn = d + f"step_{i:02d}.py"
    with open(fn) as f:
        print("="*30, fn)
        print(f.read()[:2500])
