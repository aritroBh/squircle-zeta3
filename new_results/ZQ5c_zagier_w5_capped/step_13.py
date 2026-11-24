base = '<project>/research/08_theory'
for f in ['step_09g_d13_kernel.py', 'step_12_d13_red_lll.py']:
    print("="*70)
    print(f)
    with open(f'{base}/ZQ5b_zagier_w5_larger_S/{f}') as fh:
        print(fh.read())
