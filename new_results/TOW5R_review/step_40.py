import os, time
base = '<project>/research/08_theory'
for f in ['forced_F78.pkl','tow4lib.py','result.md']:
    p = base+'/TOW4_tower_depth8/'+f
    print('TOW4', f, time.strftime('%m-%d %H:%M', time.localtime(os.path.getmtime(p))))
for f in ['tow2lib.py','L_simplified.pkl','prove_S34.py','N4_fit.json']:
    p = base+'/TOW2_tower_all_depths/'+f
    if os.path.exists(p): print('TOW2', f, time.strftime('%m-%d %H:%M', time.localtime(os.path.getmtime(p))))
# transcript: what produced N5_proved.pkl at 13:13? check transcript for 13:1
t5 = base+'/TOW5_tower_depth10'
import json
tr = open(t5+'/transcript.json').read()
i = tr.find('13:13')
print(tr[i-400:i+200] if i>0 else 'no 13:13 in transcript')
