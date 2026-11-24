base = '<project>/research/08_theory'
for f in ['proof_common.py','prove_E1_adjacent.py','prove_E2b_nonadjacent.py','prove_E2b_adjacent.py','prove_E3.py']:
    print('='*20, f)
    print(open(f'{base}/TOW1_tower_law/{f}').read())
