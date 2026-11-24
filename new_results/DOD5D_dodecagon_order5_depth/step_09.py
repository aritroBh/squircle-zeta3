base = "<project>/"
import os
for f in ["research/08_theory/A3_W4_proof/w4_symbol.py"]:
    print("="*100); print(f, os.path.getsize(base+f)); print("="*100)
    with open(base+f) as fh:
        print(fh.read())
