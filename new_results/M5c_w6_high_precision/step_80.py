import subprocess, os, time
work='<project>/research/08_theory/M5c_w6_high_precision'
time.sleep(1200)
print(subprocess.run(['bash','-lc',f'grep -E "LLL done|relations|row " "{work}/w6c_lll_x14_15_3300.out" | head -12'],capture_output=True,text=True).stdout)
pid=int(open(work+'/w6c_lll_x14_15_3300.pid').read())
try: os.kill(pid,0); print(f'pid {pid} alive')
except OSError: print(f'pid {pid} done')
print(subprocess.run(['bash','-lc','uptime'],capture_output=True,text=True).stdout)
