base = "<project>/research/08_theory"
with open(f"{base}/DCY2R_referee_varpi1/progress.log") as f:
    print("PROGRESS.LOG:"); print(f.read())
print("="*100)
for fn in ["step_49.py","step_57.py","step_58.py"]:
    print(fn); print("-"*100)
    with open(f"{base}/DCY2R_referee_varpi1/{fn}") as f:
        print(f.read())
