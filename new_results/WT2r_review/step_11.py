base = '<project>/research/08_theory/WT2b_stabilization'
for name in ['genintegrand.py', 'buildintegrand.py']:
    with open(f'{base}/{name}') as f:
        print(f"===== {name} =====")
        print(f.read())
        print()
