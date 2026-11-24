import os, json
print(open("../SRF3_octahedron_radius_proof/progress.log").read())
p = "<project>/research/10_frontier/H1_octahedron_end/direct_S_dps40.json"
print(os.path.exists(p))
if os.path.exists(p):
    d = json.load(open(p))
    for k in list(d)[:10]:
        print(k, str(d[k])[:200])
