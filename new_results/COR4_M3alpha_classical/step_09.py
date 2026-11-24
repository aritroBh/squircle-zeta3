cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
import os
for f in sorted(os.listdir(cor)):
    p = cor+f
    if f.endswith('.py'):
        with open(p) as fh:
            txt = fh.read()
        print("="*30, f, len(txt))
        print(txt[:1500])
