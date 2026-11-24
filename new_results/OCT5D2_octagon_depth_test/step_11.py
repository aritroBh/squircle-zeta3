base = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
for n in [56,57,58,59,60,61,62,63,64,65,66,67]:
    import os
    f = base+f"step_{n:02d}.py"
    if os.path.exists(f):
        print("="*30, f"step_{n:02d}.py", "="*30)
        print(open(f).read())
