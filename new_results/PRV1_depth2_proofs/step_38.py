wd = 'research/08_theory/PRV1_depth2_proofs'
for n in [29,30,31]:
    with open(f'{wd}/step_{n:02d}.py') as f:
        print('='*30, f'step_{n:02d}.py', '='*30)
        print(f.read())
