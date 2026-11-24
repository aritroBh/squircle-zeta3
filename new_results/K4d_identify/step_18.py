import os
THEORY = 'research/08_theory/'
print("08_theory subdirs:")
for d in sorted(os.listdir(THEORY)):
    p = os.path.join(THEORY,d)
    if os.path.isdir(p):
        print("  ", d)
# look for K1 folders
print("\nLooking for K1* folders anywhere in research:")
RES = 'research/'
for root,dirs,files in os.walk(RES):
    for d in dirs:
        if d.lower().startswith('k1'):
            print("  ", os.path.join(root,d))
