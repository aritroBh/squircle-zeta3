import json
j120 = json.load(open("../SRF3_octahedron_radius_proof/SRF3_hi_coeffs_K120.json"))
print(j120.keys(), j120['n'], j120['K'])
print("eta[0..5]:", j120['eta'][:6])
print("d[0..5]:", j120['d'][:6])
print("len eta:", len(j120['eta']))
j60 = json.load(open("../SRF1_radii_3d/SRF1_n3_K60.json"))
print(j60.keys())
for k in j60:
    v = j60[k]
    if isinstance(v, list):
        print(k, "len", len(v), v[:4])
    else:
        print(k, v)
