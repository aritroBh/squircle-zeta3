base = "<project>/research/08_theory/DEC1_decagon"
for s in ["step_23.py", "step_30.py", "step_31.py"]:
    print("="*30, s, "="*30)
    with open(f"{base}/{s}") as f:
        print(f.read())
