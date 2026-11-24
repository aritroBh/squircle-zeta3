import subprocess, time, os, datetime
PY='<project>/.venv/bin/python'
work='<project>/research/08_theory/M5c_w6_high_precision'
r=subprocess.run(['nice','-n','10',PY,'-B',work+'/w6c_values_extra.py','14','15','13','9','7','6'],capture_output=True,text=True,cwd=work,timeout=900)
print(r.stdout[-250:]); print('STDERR',r.stderr[-150:])
with open(work+'/progress.log','a') as f:
    f.write(f'{datetime.datetime.now().isoformat()} cumulative extras prepped: symbol stages+values for {14},{14,15},{14,15,13},{14,15,13,9},{14,15,13,9,7},{14,15,13,9,7,6} (kernels 90..95, rank stays 60; all Li6 leaves rad ~1e-3800).\n')
# check x14 LLL
print(subprocess.run(['bash','-lc',f'grep -E "LLL done|relations" "{work}/w6c_lll_x14_3300.out" | head'],capture_output=True,text=True).stdout)
pid=int(open(work+'/w6c_lll_x14_3300.pid').read())
try: os.kill(pid,0); print(f'x14 LLL pid {pid} alive')
except OSError: print(f'x14 LLL pid {pid} done')
print(subprocess.run(['bash','-lc','uptime'],capture_output=True,text=True).stdout)
