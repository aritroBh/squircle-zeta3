import subprocess
print(subprocess.run(['grep','-n','pivots\|echelon\|elim','../W5C1_certificate_scoping/progress.log'], capture_output=True, text=True).stdout[:3000])