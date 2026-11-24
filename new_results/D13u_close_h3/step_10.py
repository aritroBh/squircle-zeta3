base = '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed'
lines = open(base+'/step_38.py').read().split('\n')
for a,b in [(24,28),(49,71),(180,200)]:
    print(f"----- lines {a+1}-{b} -----")
    for i in range(a,b):
        print(f"{i+1:4d}: {lines[i]}")
