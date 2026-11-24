import subprocess, os, time
work='<project>/research/08_theory/M5c_w6_high_precision'
time.sleep(180)
print(subprocess.run(['bash','-lc',f'tail -20 "{work}/w6_vals_3500.out"'],capture_output=True,text=True).stdout)
try: os.kill(int(open(work+'/w6_vals_3500.pid').read()),0); print('vals job alive')
except OSError: print('vals job done')
