import os, glob

print("CWD:", os.getcwd())
print("\n--- CWD contents ---")
for p in sorted(os.listdir('.')):
    print(p)

print("\n--- Looking for PGEN1 / PGEN1R ---")
for root in ['.', '..', '../..', '../../..']:
    for pat in ['**/PGEN1*', '**/symbol_engine*', '**/result.md']:
        hits = glob.glob(os.path.join(root, pat), recursive=True)
        for h in hits[:50]:
            print(h)
