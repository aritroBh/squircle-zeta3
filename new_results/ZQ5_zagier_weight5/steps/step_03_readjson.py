
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
