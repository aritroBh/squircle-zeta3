for f in ["step_45.py","step_57.py","step_40.py"]:
    print("="*25, f, "="*25)
    print(open(f"../SRF3_octahedron_radius_proof/{f}").read()[:3500])
