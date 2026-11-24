base = "<project>"
for s in ["step_58.py", "step_59.py"]:
    with open(base + f"/research/08_theory/SRF1_radii_3d/{s}") as f:
        print("="*30, s, "="*30)
        print(f.read())
