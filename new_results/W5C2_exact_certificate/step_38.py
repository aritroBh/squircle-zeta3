import subprocess
print(subprocess.run(['bash','-c','ps aux | grep -v grep | grep -E "step_01|python" | head -5'], capture_output=True, text=True).stdout)
print(subprocess.run(['bash','-c','tail -c 500 step_01.out'], capture_output=True, text=True).stdout)
print(subprocess.run(['bash','-c','ls -la work/ 2>/dev/null'], capture_output=True, text=True).stdout)