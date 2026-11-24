import subprocess
result = subprocess.run(['<project>/.venv/bin/python', 
                        'research/08_theory/OCT_R34_level8/step_03_full.py'],
    capture_output=True, text=True, timeout=60, cwd='research/08_theory')
print(result.stdout)
if result.stderr:
    print("STDERR:", result.stderr)