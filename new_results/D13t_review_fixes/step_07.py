p='research/08_theory/D13s_fix_recertify/result.md'
lines=open(p).read().split('\n')
for i,l in enumerate(lines[107:116], start=108):
    print(f"{i:3d}: {l}")
print('=====')
for i,l in enumerate(lines[159:220], start=160):
    print(f"{i:3d}: {l}")
