
import json
c23 = json.load(open("../HEXC2_level6_certificate/cert_w23.json"))
c4 = json.load(open("../HEXC2_level6_certificate/cert_w4.json"))
ce = json.load(open("../HEXC2_level6_certificate/cert_euler.json"))
print("cert_w23 keys:", list(c23.keys()))
for k, v in c23.items():
    print(k, list(v.keys()), "prow len", len(v["prow"]), "y len", len(v["y"]), "y sample", v["y"][:3])
print("cert_w4 keys:", list(c4.keys()))
print("prow len", len(c4["prow"]), "y len", len(c4["y"]), "y sample", c4["y"][:3])
print("cert_euler keys:", list(ce.keys()), "prow len", len(ce["prow"]), "cols sample", ce["cols"][:5])
mm = json.load(open("../HEXC2_level6_certificate/mobius_maps.json"))
print("mobius_maps.json:", str(mm)[:800])
