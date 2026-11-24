import json
out = json.load(open("<project>/research/10_frontier/E1_edge_law/verify/out/t3_b0_150.json"))
print(json.dumps(out, indent=1)[:4000])
