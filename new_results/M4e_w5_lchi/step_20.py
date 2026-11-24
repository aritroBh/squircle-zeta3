import json
p = '<project>/research/08_theory/M4d_w5_nice'
ksat = json.load(open(p + '/nice_ksat.json'))
print("names:", ksat['names'])
print("num kernel rows:", len(ksat['K']))
for row in ksat['K']:
    print(row)
