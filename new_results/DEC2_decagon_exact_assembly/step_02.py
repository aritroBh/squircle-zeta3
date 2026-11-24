import os
base = "<project>/research/08_theory"
for d in ["DEC2_decagon_exact_assembly", "DEC1_decagon", "PRV1_depth2_proofs",
          "OCT5_octagon_exact", "OCT6_octagon_closed", "CX1_octagon_depth2_proof"]:
    p = os.path.join(base, d)
    print("="*80)
    print(p, "exists:", os.path.isdir(p))
    if os.path.isdir(p):
        for f in sorted(os.listdir(p)):
            full = os.path.join(p, f)
            print(f"   {f}  ({os.path.getsize(full)} bytes)")
