base = "<project>/research/08_theory/DEC1_decagon"
for s in ["step_33.py", "step_34.py"]:
    print("="*30, s, "="*30)
    with open(f"{base}/{s}") as f:
        print(f.read())
