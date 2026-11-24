print(open('../W5C1_certificate_scoping/w5_fold2.py').read())
print("===== progress.log tail =====")
import subprocess
print(subprocess.run(['tail','-40','../W5C1_certificate_scoping/progress.log'], capture_output=True, text=True).stdout)