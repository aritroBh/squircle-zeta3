import subprocess, glob
files = sorted(glob.glob('../OCT9b_octagon_debug/step_*.py'))
for pat in ['solB0','sol_B123','W_parts','brackets.pkl','Sunit_lists','sols2357','lam_cache','reduce_e1_e2']:
    hits = [f for f in files if pat in open(f).read()]
    print(pat, '->', hits)
