
import json
s3 = '<project>/research/07_push/S3_superellipsoid/claude'
sv = json.load(open(s3+'/S_values.json'))
print("S_values.json:", type(sv), list(sv)[:10] if isinstance(sv,dict) else len(sv))
s = json.dumps(sv)[:600]; print(s)
svh = json.load(open(s3+'/S_values_huge.json'))
print("S_values_huge.json:", type(svh), list(svh)[:10] if isinstance(svh,dict) else len(svh))
print(json.dumps(svh)[:300])
v2 = json.load(open(s3+'/../verify/v02_coeffs_h64.json'))
print("v02_coeffs_h64:", type(v2), list(v2)[:6] if isinstance(v2,dict) else len(v2))
print(json.dumps(v2)[:400])
