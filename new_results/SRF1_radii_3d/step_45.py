
import json
vf = '<project>/research/07_push/S3_superellipsoid/verify'
v2 = json.load(open(vf+'/v02_coeffs_h64.json'))
for key in v2:
    val = v2[key]
    if isinstance(val, list):
        print(key, "list len", len(val), "first:", str(val[0])[:80], " last:", str(val[-1])[:80])
    else:
        print(key, val)
print()
print(open(vf+'/v02_coeffs.py').read()[:3500])
