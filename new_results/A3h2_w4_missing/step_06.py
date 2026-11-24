base='<project>'
import subprocess
print(subprocess.run(['ls','-la',base+'/research/08_theory/A3h_weight4_target/Q'],capture_output=True,text=True).stdout)
print(subprocess.run(['ls','-la',base+'/research/08_theory/A3c_claude_elementary/Q'],capture_output=True,text=True).stdout)
