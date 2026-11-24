import os, sys, subprocess
src = r'''
import os, json
d = "../M4i_w5_noL"
rm = open(os.path.join(d,"result.md")).read()
print("result.md length:", len(rm))
import re
# print planner-ish and verification sections
idx = rm.lower().find("planner")
print("---- around 'planner' ----")
print(rm[max(0,idx-500):idx+2500] if idx>=0 else "(no planner mention)")
print("---- head ----")
print(rm[:1500])
for fn in ["Lchi8_5.json"]:
    J = json.load(open(os.path.join(d,fn)))
    s = json.dumps(J)
    print(fn, ":", s[:600])
'''
open('steps/step_04_readmd.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_04_readmd.py'],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-800:])
with open('progress.log','a') as f: f.write("step_04_readmd: read M4i result.md\n")
