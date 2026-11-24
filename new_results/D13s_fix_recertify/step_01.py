
import os, shutil, subprocess

BASE = 'research/08_theory'
WORK = os.path.join(os.getcwd(), 'D13s_work')
os.makedirs(WORK, exist_ok=True)
print("WORK =", WORK)

# Copy the relevant directories/files into the work directory
srcs = [
    'D13g_all_N/step_38.py',
    'D13h_final_cell/step_10_J192.py',
    'D13h_final_cell/step_18.py',
    'D13h_final_cell/planner_sweep_ratio.py',
]
for s in srcs:
    src = os.path.join(BASE, s)
    dst = os.path.join(WORK, s.replace('/', '__'))
    shutil.copy(src, dst)
    print("copied", s, "->", dst, os.path.getsize(dst), "bytes")

# Also copy the review result.md for reading
rev = os.path.join(BASE, 'D13r_lemma_review/result.md')
shutil.copy(rev, os.path.join(WORK, 'review__result.md'))
print("copied review result.md", os.path.getsize(rev), "bytes")

# also check for review progress.log
rlog = os.path.join(BASE, 'D13r_lemma_review/progress.log')
if os.path.exists(rlog):
    shutil.copy(rlog, os.path.join(WORK, 'review__progress.log'))
    print("copied review progress.log", os.path.getsize(rlog), "bytes")
