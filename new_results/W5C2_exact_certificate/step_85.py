import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 03:48:57 step_01b pp DONE: |S|=12156=rank (engine 12120 + framew5 36, NO frame234), quotient 62; skip-ahead at row 27820; D_rat & D_sq in span mod p1; 679s.\n')
print(subprocess.run(['bash','-c','grep -E "pmb|step_01b pm" step_01b.out | tail -6'], capture_output=True, text=True).stdout)
time.sleep(420)
print(subprocess.run(['bash','-c','grep -E "pmb|step_01b pm|mpb|step_01b mp" step_01b.out | tail -10'], capture_output=True, text=True).stdout)