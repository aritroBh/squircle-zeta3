import subprocess, time, os
work='<project>/research/08_theory/M5c_w6_high_precision'
time.sleep(300)
r=subprocess.run(['bash','-lc',f'grep -E "bits=11800|complete|agreement|absdiff" "{work}/w6_hp.out" | cut -c1-200'],capture_output=True,text=True)
print(r.stdout)
try: os.kill(6401,0); print('pid 6401 alive')
except OSError: print('pid 6401 done')
