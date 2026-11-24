base = "<project>"
for s in ["step_55.py", "step_56.py", "step_57.py"]:
    with open(base + f"/research/08_theory/SRF1_radii_3d/{s}") as f:
        print("="*30, s, "="*30)
        print(f.read())
