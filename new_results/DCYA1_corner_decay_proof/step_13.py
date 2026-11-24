import subprocess, os
for d in ['D13d_decay_theorem','D13s_fix_recertify','D13u_close_h3','K2_corner_decay']:
    base=f'<project>/research/08_theory/{d}/'
    print(f"\n\n############ {d} ############")
    print(subprocess.run(['ls',base], capture_output=True, text=True).stdout)
