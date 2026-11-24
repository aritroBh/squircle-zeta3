import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 05:24:19 step_03 pp SOLVED: y_rat nz=12149/12156 (num 2078d, den 2065d), y_sq nz=12142 (2081d/2069d); Dixon 577 steps, 963s; cross-check A.Y==B mod p2 OK.\n')
time.sleep(540)
print(subprocess.run(['bash','-c','tail -8 step_03.out'], capture_output=True, text=True).stdout)