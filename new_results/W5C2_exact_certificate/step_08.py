import subprocess
print(subprocess.run(['ls','../W5C1_certificate_scoping/sector/'], capture_output=True, text=True).stdout)