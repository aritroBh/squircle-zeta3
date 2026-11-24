import subprocess, os
os.chdir('<project>')
print(os.getcwd())
print(subprocess.run(['ls','-la'],capture_output=True,text=True).stdout)
print("--- CHN4 dir ---")
print(subprocess.run(['ls','-la','CHN4_p5_w11_locus'],capture_output=True,text=True).stdout)
