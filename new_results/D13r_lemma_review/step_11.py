base='research/08_theory'
src=open(base+'/D13h_final_cell/step_10_J192.py').read()
i=src.find('cfp')
print(src[i:i+3300])
