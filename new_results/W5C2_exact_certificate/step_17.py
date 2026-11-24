import subprocess
out = subprocess.run(['grep','-n','sector\|membership\|in span\|D_rat\|D_sq','../W5C1_certificate_scoping/progress.log'], capture_output=True, text=True).stdout
print(out[:4000])