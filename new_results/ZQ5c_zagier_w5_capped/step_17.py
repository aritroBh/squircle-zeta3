base = '<project>/research/08_theory'
for f in ['step_05_d17_final.py', 'step_10_d13_vals.py']:
    print("="*70)
    print(f)
    with open(f'{base}/ZQ5b_zagier_w5_larger_S/{f}') as fh:
        print(fh.read())
