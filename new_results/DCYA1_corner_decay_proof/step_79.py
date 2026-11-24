import subprocess, os
for d in ['D13s_fix_recertify/D13s_work','D13u_close_h3','D13d_decay_theorem']:
    base=f'<project>/research/08_theory/{d}'
    print(f"=== {d} ===")
    print(subprocess.run(['ls',base],capture_output=True,text=True).stdout)
    # look for fixed/ subdirectory
    import glob
    for sub in glob.glob(base+'/**/step_38*.py',recursive=True)+glob.glob(base+'/**/planner_sweep*.py',recursive=True):
        print("  found:",sub)
