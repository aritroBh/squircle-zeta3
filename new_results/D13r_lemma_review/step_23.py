base='research/08_theory'
for f in ['step_08.py','step_09.py']:
    print('='*20, f, '='*20)
    print(open(base+'/D13h_final_cell/'+f).read()[:4200])
