
import subprocess
print(subprocess.run(['find','<project>','-name','main.tex','-not','-path','*/node_modules/*'],capture_output=True,text=True).stdout)
