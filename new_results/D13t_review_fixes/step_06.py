p='research/08_theory/D13s_fix_recertify/result.md'
lines=open(p).read().split('\n')
print(len(lines))
for i,l in enumerate(lines[55:160], start=56):
    print(f"{i:3d}: {l}")
