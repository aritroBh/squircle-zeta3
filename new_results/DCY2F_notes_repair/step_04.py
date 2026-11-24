import os
theory = "<project>/research/08_theory"
# Find DCY-related dirs
for item in sorted(os.listdir(theory)):
    if "DCY" in item.upper():
        p = os.path.join(theory, item)
        print(f"  {'DIR ' if os.path.isdir(p) else 'FILE'} {item}")
