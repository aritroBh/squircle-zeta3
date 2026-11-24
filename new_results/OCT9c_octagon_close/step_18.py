import subprocess
for pat in ['solB0','sol_B','W_parts','brackets','Sunit_lists','sols2357','lam_cache']:
    out = subprocess.run(['grep','-l',pat,'../OCT9b_octagon_debug/step_7'+str(i)+'.py' for i in range(10)],capture_output=True,text=True)
    print(pat, '->', out.stdout.split())
