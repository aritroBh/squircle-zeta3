import os
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
summary = """# W5C2 — exact certificate for W5 = N5 + R5_min3 (status: RUNNING, step 1 nearly done)

Pipeline (mirrors A3h5 at weight 5):
- step_01b: row trimming per sector (LeanEchelon sparse GE mod p1=16777213, plateau skip-ahead).
  DONE pp: |S|=12156=rank (engine 12120 + framew5 36), quotient 62, 679s.
  DONE pm: |S|=11101=rank (engine 11065 + framew5 36), quotient 63, 992s.
  RUNNING mp, mm. NOTE: D needs NO frame234 rows in pp/pm (engine+36 framew5 suffice).
- step_02: second-prime (p2=16777199) checks: rank(S)=|S|, det(R_{S,C})!=0, D in span(S).
- step_03: exact solve y.R_{S,C}=D_C over Q per sector (custom Dixon, both RHS).
- step_04: exact Fraction verification on all sector columns + fold/unfold lemma; certificates.
- verify_W5.py: standalone re-check from saved files.
"""
open('summary.md','w').write(summary)
print('summary updated')
import subprocess, time
time.sleep(500)
print(subprocess.run(['bash','-c','grep -E "mpb|step_01b mp|mmb|step_01b mm|plateau|feed done" step_01b.out | tail -12'], capture_output=True, text=True).stdout)