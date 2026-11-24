base = '<project>/research/08_theory/'
pg = base+'PGEN1R_referee_alphabet/'
import os
# grep for delta in the step files
for f in sorted(os.listdir(pg)):
    if f.startswith('step_') and f.endswith('.py'):
        c = open(pg+f).read()
        if 'delta' in c or 'wedge' in c or 'N3' in c or 'Xi' in c or '(2,1)' in c:
            print(f"--- {f} ---")
            print(c[:1600])
            print()
