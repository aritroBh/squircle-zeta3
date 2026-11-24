import subprocess
print(open('../W5C1_certificate_scoping/step_113.py').read())
print(subprocess.run(['bash','-c','ps aux | grep -v grep | grep step_01'], capture_output=True, text=True).stdout)