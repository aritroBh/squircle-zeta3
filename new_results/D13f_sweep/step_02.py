
import os
os.chdir("<home>)")
# find prev_result.md, step_34.py, step_36.py, mk_deg9.tsv anywhere
import subprocess
for target in ["prev_result.md", "step_34.py", "step_36.py", "mk_deg9.tsv"]:
    for root, dirs, files in os.walk("."):
        # skip .venv
        if ".venv" in root:
            continue
        for f in files:
            if f == target:
                print(os.path.join(root, f))
