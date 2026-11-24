import subprocess
result = subprocess.run(['<project>/.venv/bin/python', 
                        'research/08_theory/OCT_R34_level8/step_04_lam_proof.py'],
    capture_output=True, text=True, timeout=120, cwd='research/08_theory')
print(result.stdout)
if result.stderr:
    print("\nSTDERR:", result.stderr)