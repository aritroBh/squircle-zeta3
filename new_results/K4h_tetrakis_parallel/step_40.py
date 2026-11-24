import json
r = json.load(open("ckpt_h768/chunk_0000.json"))
print("chunk_0000 keys:", list(r.keys()))
print("idx:", r['idx'], "n:", r['n'], "el:", round(r['el'],1), "stats:", r['stats'])
print("acc digits:", [len(a.split('.')[-1].rstrip('eE0123456789') if False else a) for a in r['acc']][:1])
print("acc[0] head:", r['acc'][0][:60], "... len:", len(r['acc'][0]))
r1 = json.load(open("ckpt_h768/chunk_0001.json"))
print("chunk_0001 stats:", r1['stats'], "el:", round(r1['el'],1))
