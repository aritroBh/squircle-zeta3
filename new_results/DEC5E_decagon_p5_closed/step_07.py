import subprocess
# look at the step_50..52 files in DEC5D
for i in [50,51,52]:
    print('='*90)
    print(f'step_{i}.py')
    print(open(f'../DEC5D_decagon_order5_depth/step_{i:02d}.py').read())
