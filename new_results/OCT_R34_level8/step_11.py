import subprocess
import os
# Run from the 08_theory directory needed for imports
result = subprocess.run(['<project>/.venv/bin/python', 
                        'research/08_theory/OCT_R34_level8/step_01_atoms.py'],
    capture_output=True, text=True, timeout=60, cwd='research/08_theory')
print(result.stdout)
if result.stderr:
    print("STDERR:", result.stderr)