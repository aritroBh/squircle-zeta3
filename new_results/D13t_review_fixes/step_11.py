p='research/08_theory/D13s_fix_recertify/D13s_work/fixed/step_10_J192.py'
lines=open(p).read().split('\n')
for i,l in enumerate(lines[:150], start=1):
    print(f"{i:3d}: {l}")
