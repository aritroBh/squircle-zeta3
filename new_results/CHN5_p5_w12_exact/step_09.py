base='<project>/research/08_theory'
import subprocess
print(subprocess.run(['ls','-la',base+'/CHN4_p5_w11_locus/planner'],capture_output=True,text=True).stdout)
print(subprocess.run(['sed','-n','1,120p',base+'/CHN4_p5_w11_locus/planner/planner_check.out'],capture_output=True,text=True).stdout if os.path.exists(base+'/CHN4_p5_w11_locus/planner/planner_check.out') else "no planner_check.out")
import os
for root,d,fs in os.walk(base+'/CHN4_p5_w11_locus/planner'):
    for f in fs:
        print(os.path.join(root,f))
