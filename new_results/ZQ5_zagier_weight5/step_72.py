import sys, subprocess
src = r'''
import sys; sys.path.insert(0,'.')
import w5lib as W
import pipe
C = W.letters_info(13)
for B in (10, 14):
    cands = W.enum_exceptional(C, B=B)
    print(f"B={B}: {len(cands)} exceptional S-units")
    if len(cands) > 4:
        for c in cands:
            print(f"   x={pipe.fmt_qd(c['x']):>22s} vx={c['vx']} 1-x={c['v1m']} h={c['h']}")
'''
open('steps/step_45_d13_more.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_45_d13_more.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1000:])
with open('progress.log','a') as f: f.write("step_45_d13_more: larger B for D=13\n")
