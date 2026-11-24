base = '<project>'
for f in ['progress.log', 'step_03.py', 'step_08.py']:
    print('='*100); print(f); print('='*100)
    with open(base + '/research/08_theory/ZLAD1_zeta_ladder/' + f) as fh:
        print(fh.read())
