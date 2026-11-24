# Copy the script to the correct directory
import shutil
shutil.copy('step_01_atoms.py', 'research/08_theory/step_01_atoms.py')
result = subprocess.run(['<project>/.venv/bin/python', 'step_01_atoms.py'],
    capture_output=True, text=True, timeout=60, cwd='research/08_theory')
print(result.stdout)
if result.stderr:
    print("STDERR:", result.stderr)