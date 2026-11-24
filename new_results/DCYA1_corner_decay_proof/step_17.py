import subprocess, os
for d in ['K3_fourier_all_orders','AH3_fn_any_angle','D13c_two_ends']:
    base=f'<project>/research/08_theory/{d}/'
    rm=base+'result.md'
    if os.path.exists(rm):
        print(f"\n\n######## {d}/result.md (first 3500) ########")
        print(subprocess.run(['head','-c','3500',rm], capture_output=True, text=True).stdout)
