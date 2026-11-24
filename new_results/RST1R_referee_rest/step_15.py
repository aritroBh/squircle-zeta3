import subprocess
for d in ['DCY2_decay_first_correction','CX6R_review','DCY2R_referee_varpi1','RST1_rest_exponent/planner']:
    p=f'<project>/research/08_theory/{d}'
    print("="*20, d)
    print(subprocess.run(['ls','-la',p], capture_output=True, text=True).stdout)
