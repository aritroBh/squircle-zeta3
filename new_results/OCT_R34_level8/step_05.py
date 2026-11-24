# Now let's run the template to see it work (from the correct directory)
import subprocess
result = subprocess.run([
    '<project>/.venv/bin/python',
    'research/08_theory/SQ_chi3_bridge/verify_bridge_claude.py'
], cwd='research/08_theory',
    capture_output=True, text=True, timeout=120)
print("STDOUT:")
print(result.stdout)
print("\nSTDERR:")
print(result.stderr)