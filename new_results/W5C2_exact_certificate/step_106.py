import os, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 05:08:00 BUGFIX: step_03/step_04 loaded stale work/rows_*.pkl (v1 ordering) instead of work/feed_*.pkl (step_01b ordering that S indexes into) -> singular A. Fixed; deleted stale rows_*.pp solve restarted. nmod inv(12156)=406s measured; fmpz solve rejected (148s at n=4000).\n')
time.sleep(450)
print(subprocess.run(['bash','-c','tail -8 step_03.out'], capture_output=True, text=True).stdout)