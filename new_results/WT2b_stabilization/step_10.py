src='<project>/research/08_theory/WT2_weight2_all_k'
for f in ['step_48.py','step_49.py','step_50.py','step_51.py','step_52.py','step_53.py','step_54.py']:
    c = open(f'{src}/{f}').read()
    print('='*20, f, len(c))
    print(c[:1800])
