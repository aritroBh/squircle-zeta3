import subprocess, time, os
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
r=subprocess.run(['nice','-n','10',PY,'-B',work+'/w6c_values_extra.py','14','15','13'],capture_output=True,text=True,cwd=work,timeout=900)
print(r.stdout[-300:]); print('STDERR',r.stderr[-200:])
# check the {14} LLL
print(subprocess.run(['bash','-lc',f'grep -E "LLL done|relations" "{work}/w6c_lll_x14_3300.out" | head'],capture_output=True,text=True).stdout)
pid=int(open(work+'/w6c_lll_x14_3300.pid').read())
try: os.kill(pid,0); print(f'x14 LLL pid {pid} alive')
except OSError: print(f'x14 LLL pid {pid} done')
