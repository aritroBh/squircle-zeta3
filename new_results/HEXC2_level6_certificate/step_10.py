import os, shutil
base = '<project>/research/08_theory'
cwd = os.getcwd()
# symlink-free: we'll just sys.path to inputs. First check the paper theorem reference
import subprocess
print(subprocess.run(['grep','-n','thm:W4R','-r','<project>/paper/ieee/main.tex'],capture_output=True,text=True).stdout[:1500])
