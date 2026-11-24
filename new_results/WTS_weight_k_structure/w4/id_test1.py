import sys, json
from id_lib import *
dps = 300; set_dps(dps); print("mp.mp.dps =", mp.mp.dps)
I = load_I(dps)
c = [combo(I, f"c{j}") for j in range(3)]
ref = json.load(open(os.path.join(RES, "id_coeffs_dps300.json")))["c"]
for j in range(3): print(f"c{j}(from I_n) - c{j}(fit) =", mp.nstr(c[j] - mp.mpf(ref[j]), 5))
# k=3 check
print("W_4(3) from I_n - fit:", mp.nstr(combo(I, 'k3') - mp.mpf(json.load(open(os.path.join(RES,'id_coeffs_dps300.json')))['W4']['3']), 5))
