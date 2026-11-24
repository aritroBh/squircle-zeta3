import json
b3 = "<project>/research/07_push/S3_superellipsoid/claude/"
fc = json.load(open(b3+"final_constants_400.json"))
for k in fc: 
    v = fc[k]
    print(k, (v[:60] if isinstance(v,str) else v))
