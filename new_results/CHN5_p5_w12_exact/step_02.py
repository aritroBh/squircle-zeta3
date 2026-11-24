import subprocess, os
os.chdir('<project>')
print(subprocess.run(['find','.','-maxdepth','3','-name','*CHN4*','-o','-maxdepth','3','-name','*p5*'],capture_output=True,text=True).stdout)
print(subprocess.run(['ls','research'],capture_output=True,text=True).stdout)
