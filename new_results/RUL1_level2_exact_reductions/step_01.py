import subprocess
print(subprocess.run(['bash','-c','pwd; ls -la; echo ---; cat progress.log 2>/dev/null'], capture_output=True, text=True).stdout)
