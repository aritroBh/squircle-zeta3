base = "<project>/research/08_theory/OCT6_octagon_closed"
for s in ["step_17.py", "step_18.py"]:
    print("="*30, s, "="*30)
    with open(f"{base}/{s}") as f:
        print(f.read())
