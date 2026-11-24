
base = '<project>/research/08_theory'
for f in ['K1f_hexagon_order5/step_01.py', 'K1f_hexagon_order5/step_02.py']:
    print(f'===== {f} =====')
    with open(f'{base}/{f}') as fh:
        print(fh.read())
    print()
