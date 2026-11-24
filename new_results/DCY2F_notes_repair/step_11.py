base = "<project>/research/08_theory"
for fn in ["step_06_finite.py", "step_02_varpi1.py", "step_03_watson.py", "step_07_partb.py"]:
    print("="*100)
    print(fn)
    print("="*100)
    with open(f"{base}/DCY2_decay_first_correction/{fn}") as f:
        print(f.read())
