base = "<project>/research/08_theory"
for fn in ["dcy2_sweep.py", "dcy2_regionc.py", "dcy2_common.py"]:
    print("="*100); print(fn); print("="*100)
    with open(f"{base}/DCY2_decay_first_correction/{fn}") as f:
        print(f.read())
