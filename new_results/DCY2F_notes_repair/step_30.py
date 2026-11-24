base = "<project>/research/08_theory"
with open(f"{base}/DCY2_decay_first_correction/step_02_varpi1.py") as f:
    print(f.read())
print("#"*100)
with open(f"{base}/DCY2_decay_first_correction/dcy2_audit.py") as f:
    print(f.read()[:1500])
