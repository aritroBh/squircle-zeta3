import datetime, subprocess, os
work='<project>/research/08_theory/M5c_w6_high_precision'
def log(s):
    with open(work+'/progress.log','a') as f:
        f.write(f'{datetime.datetime.now().isoformat()} {s}\n')
log('M5c start: copied M5b scripts+data; flint planted LLL control PASS (gram=exact)')
log('Lchi8_5 at 3600 digits (flint Hurwitz, 12600 bits, rad 5.3e-3793); vs Lchi8_5_2500 absdiff 3.4e-2511; saved Lchi8_5_3600.json')
log('launched W6 quadrature 11500+11800 bits (m5b_precision.py w6), pid 6401')
# check job
r=subprocess.run(['bash','-lc',f'cat "{work}/w6_hp.out"; echo ===; grep -c . "{work}/progress.log"'],capture_output=True,text=True)
print(r.stdout)
try: os.kill(6401,0); print('pid 6401 alive')
except OSError: print('pid 6401 done')
