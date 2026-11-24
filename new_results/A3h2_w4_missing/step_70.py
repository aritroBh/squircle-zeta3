import subprocess
out = subprocess.run(['find','<project>','-name','qsqrt2*','-not','-path','*/node_modules/*'],capture_output=True,text=True)
print(out.stdout, out.stderr)
