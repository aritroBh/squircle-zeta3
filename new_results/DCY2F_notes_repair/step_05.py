import os
for d in ["DCY2_decay_first_correction", "DCY2R_referee_varpi1"]:
    root = os.path.join("<project>/research/08_theory", d)
    print("="*80)
    print(root)
    for dirpath, dirnames, filenames in os.walk(root):
        for f in sorted(filenames):
            p = os.path.join(dirpath, f)
            print(f"  {os.path.relpath(p, root)}  ({os.path.getsize(p)} bytes)")
