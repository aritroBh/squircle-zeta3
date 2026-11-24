base = "<project>/research/08_theory/DEC1_decagon"
for s in ["step_36.py", "step_37.py"]:
    print("="*30, s, "="*30)
    with open(f"{base}/{s}") as f:
        print(f.read())
