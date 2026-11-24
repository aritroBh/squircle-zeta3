import os
os.chdir('research/08_theory/DEC1_decagon')
files = sorted(f for f in os.listdir('.') if not f.startswith('.'))
with open('progress.log','a') as f:
    f.write("step: all four orbit equations verified 830 dps. Work complete: P'''(1)=60*L*S (15 terms),\n"
            " verified 1.1e-676; no depth-2 constant survives; local z3 == C_5 exactly. result.md final.\n")
print("files in DEC1_decagon:")
for f in files: print(" ", f)