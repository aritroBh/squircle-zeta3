import subprocess, os, sys

print("Running verification script...")
result = subprocess.run([sys.executable, 'verify_kspan.py'], 
                       capture_output=True, text=True)
print("STDOUT:")
print(result.stdout)
if result.stderr:
    print("STDERR:")
    print(result.stderr)
print(f"Return code: {result.returncode}")