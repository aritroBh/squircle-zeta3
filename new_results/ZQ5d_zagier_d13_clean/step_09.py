import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
gl = json.load(open(c+'/out/d13/genuine_lattice_cap.json'))
print(type(gl), list(gl.keys()) if isinstance(gl,dict) else len(gl))
print(json.dumps(gl, indent=1)[:4000])
