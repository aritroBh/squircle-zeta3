import subprocess
print(subprocess.run(['bash','-c','cat step_01.py; echo ===; ls <project>/research/08_theory/ | head -80'], capture_output=True, text=True).stdout)
