import subprocess
print(subprocess.run(['grep','-rln','rows_framew5','../W5C1_certificate_scoping/'], capture_output=True, text=True).stdout)
print(subprocess.run(['ps','aux'], capture_output=True, text=True).stdout[:200])
print(subprocess.run(['bash','-c','ps aux | grep -v grep | grep step_01'], capture_output=True, text=True).stdout)