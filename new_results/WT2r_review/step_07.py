base = '<project>/research/08_theory/WT2b_stabilization'
for n in [40, 41]:
    with open(f'{base}/step_{n}.py') as f:
        print(f"===== step_{n}.py =====")
        print(f.read())
        print()
