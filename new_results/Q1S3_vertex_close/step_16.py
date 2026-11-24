import os, json
print("cwd:", os.getcwd())
print(os.listdir('.'))
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
r1 = json.load(open(b3+"V1_1250.json"))
r3 = json.load(open(b3+"V3p_1250.json"))
print(r1.keys(), r3.keys())
print(r1['V1'][:80])
print(r3['V3p'][:80])
# check other json files for Q1 pieces values
print(json.load(open(b3+"Q1_sub_400.json")).keys() if os.path.exists(b3+"Q1_sub_400.json") else "no Q1_sub_400")
