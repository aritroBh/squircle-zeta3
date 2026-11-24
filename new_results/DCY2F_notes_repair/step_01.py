import os
# Explore cwd and readable dirs
for root in ["research/08_theory/DCY2F_notes_repair",
             "research/08_theory/DCY2_decay_first_correction",
             "research/08_theory/DCY2R_referee_varpi1"]:
    print("="*80)
    print(root, "EXISTS:", os.path.exists(root))
    if os.path.exists(root):
        for dirpath, dirnames, filenames in os.walk(root):
            for f in sorted(filenames):
                p = os.path.join(dirpath, f)
                print(f"  {p}  ({os.path.getsize(p)} bytes)")
