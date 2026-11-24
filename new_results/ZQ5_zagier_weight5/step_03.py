import os, sys, subprocess
src = r'''
import os, json
d = "../M4i_w5_noL"
print("M4i dir listing:", sorted(os.listdir(d)))
J = json.load(open(os.path.join(d,"L_relation.json")))
print("json type:", type(J))
if isinstance(J, dict):
    print("keys:", list(J.keys()))
s = json.dumps(J)
print("json length:", len(s))
print(s[:3000])
'''
open('steps/step_03_readjson.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_03_readjson.py'],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-800:])
with open('progress.log','a') as f: f.write("step_03_readjson: inspected L_relation.json\n")
