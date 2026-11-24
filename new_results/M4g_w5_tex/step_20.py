import subprocess
out=subprocess.run(['grep','-rn','-e','cobracket','-e','Lie word','-e','\\[X','../M4d_w5_nice/result.md','../M4c_w5_closed_form/result.md'],capture_output=True,text=True).stdout
print(out[:4000])
