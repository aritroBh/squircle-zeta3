import subprocess, time
time.sleep(120)
print(subprocess.run(['bash','-c','grep -E "step_02 mm|DONE" step_02.out | tail -6'], capture_output=True, text=True).stdout)
# log progress
with open('<project>/research/08_theory/W5C2_exact_certificate/progress.log','a') as f:
    f.write('2026-10-29 04:37:00 step_02: pp/pm/mp ALL CHECKS PASSED mod p2 (rank |S|, D in span, det(R_S,C)!=0). mm pending.\n')