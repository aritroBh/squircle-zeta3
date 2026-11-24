p='research/08_theory/D13s_fix_recertify/D13s_work/fixed/step_38.py'
lines=open(p).read().split('\n')
for i,l in enumerate(lines[78:152], start=79):
    print(f"{i:3d}: {l}")
