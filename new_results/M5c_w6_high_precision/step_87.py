import subprocess, os, time
work='<project>/research/08_theory/M5c_w6_high_precision'
time.sleep(1080)
print(subprocess.run(['bash','-lc',f'cat "{work}/w6c_control.out"'],capture_output=True,text=True).stdout)
pid=int(open(work+'/w6c_control.pid').read())
try: os.kill(pid,0); print(f'pid {pid} alive')
except OSError: print(f'pid {pid} done')
print(subprocess.run(['bash','-lc','uptime'],capture_output=True,text=True).stdout)
