import subprocess
print(subprocess.run(['grep','-ln','rows_framew5','../W5C1_certificate_scoping/'], capture_output=True, text=True).stdout)
print(subprocess.run(['grep','-rn','rows_framew5','../W5C1_certificate_scoping/step_11*.py'], capture_output=True, text=True).stdout[:1500])
print(subprocess.run(['tail','-4','step_01.out'], capture_output=True, text=True).stdout)