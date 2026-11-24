import subprocess, os, time
work='<project>/research/08_theory/M5c_w6_high_precision'
time.sleep(480)
print(subprocess.run(['bash','-lc',f'cat "{work}/w6c_lll_3000.out"; echo ===; tail -4 "{work}/progress.log"'],capture_output=True,text=True).stdout)
pid=int(open(work+'/w6c_lll_3000.pid').read())
try: os.kill(pid,0); print(f'pid {pid} alive')
except OSError: print(f'pid {pid} done')
