import subprocess, os, time
work='<project>/research/08_theory/M5c_w6_high_precision'
time.sleep(1500)
print(subprocess.run(['bash','-lc',f'grep -E "LLL done|relations" "{work}/w6c_lll_x14_15_9_3300.out" | head'],capture_output=True,text=True).stdout)
pid=int(open(work+'/w6c_lll_x14_15_9_3300.pid').read())
try: os.kill(pid,0); print(f'pid {pid} alive')
except OSError: print(f'pid {pid} done')
print(subprocess.run(['bash','-lc','uptime'],capture_output=True,text=True).stdout)
