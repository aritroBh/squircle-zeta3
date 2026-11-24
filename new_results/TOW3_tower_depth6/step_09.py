base = '<project>/research/08_theory'
work = f'{base}/TOW3_tower_depth6'
with open(f'{work}/progress.log','a') as f:
    f.write('step_01-07: read TOW2 result/proof/tow2lib; copied library+pickles+proof scripts into TOW3\n')
print(open(f'{base}/TOW2_tower_all_depths/prove_S34.py').read()[:3500])
