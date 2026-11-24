import os
# Look at the parent theory directory
theory = "<project>/research/08_theory"
print("08_theory contents:")
for item in sorted(os.listdir(theory)):
    p = os.path.join(theory, item)
    print(f"  {'DIR ' if os.path.isdir(p) else 'FILE'} {item}")
