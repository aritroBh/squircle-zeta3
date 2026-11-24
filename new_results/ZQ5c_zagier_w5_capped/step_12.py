base = '<project>/research/08_theory'
for f in ['step_06_w3_d5.py', 'step_07_w3_d13.py']:
    print("="*70)
    print(f)
    with open(f'{base}/ZQ5b_zagier_w5_larger_S/{f}') as fh:
        print(fh.read())
