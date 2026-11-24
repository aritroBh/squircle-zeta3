
import os, json
base = "<project>"
d = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/direct_S_dps40.json")))
print(json.dumps(d, indent=1)[:2500])
