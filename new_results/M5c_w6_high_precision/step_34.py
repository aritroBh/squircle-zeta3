import json, datetime
import mpmath as mp
work='<project>/research/08_theory/M5c_w6_high_precision'
mp.mp.dps=1800
d12=json.load(open(work+'/w6_11800.json'))
w54=json.load(open(work+'/W6_5400.json'))
diff=abs(mp.mpf(d12['value_full'])-mp.mpf(w54['value']))
print('W6(11800 bits) vs W6_5400 absdiff:', mp.nstr(diff,8))
def log(s):
    with open(work+'/progress.log','a') as f:
        f.write(f'{datetime.datetime.now().isoformat()} {s}\n')
log('W6 11500 bits done (262s): radius 5.0e-3459, tail 1.0e-3474, 3431 digits')
log('W6 11800 bits done (553s): radius 2.6e-3549, tail 5.7e-3565, 3522 digits; agreement 11500/11800 = 2.3e-3432; vs W6_400 3.7e-401')
log(f'control W6(11800) vs W6_5400 absdiff = {mp.nstr(diff,8)}')
print('logged')
