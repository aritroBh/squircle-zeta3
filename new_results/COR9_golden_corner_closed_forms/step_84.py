import subprocess
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
proc = subprocess.run(['<project>/.venv/bin/python', WD+'verify_cor9.py'],
                      capture_output=True, text=True, timeout=280)
for line in proc.stdout.splitlines():
    if 'FAIL' in line: print(line)
