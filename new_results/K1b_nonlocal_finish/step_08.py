import subprocess
d = "<project>/research/08_theory/K1_nonlocal_law/"
for i in range(28, 34):
    fn = d + f"step_{i:02d}.py"
    try:
        with open(fn) as f:
            head = f.read()
        print("="*30, fn)
        print(head[:600])
    except FileNotFoundError:
        print(fn, "missing")
