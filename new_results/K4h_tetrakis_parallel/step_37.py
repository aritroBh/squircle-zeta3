import json
import mpmath as mp
mp.mp.dps = 100
h96 = json.load(open("<project>/research/08_theory/K4e_tetrakis_d3/modelA_h96.json"))
val = json.load(open("val96_result.json"))
print("val96 parallel vs K4e modelA_h96.json (same h=1/96 grid, independent WMIN/dps):")
for k in range(5):
    a = mp.mpf(val["Spar"][k]); b = mp.mpf(h96["Sdk"][k])
    print("  k=%d  |diff| = %s" % (k, mp.nstr(abs(a-b),3)))
print("stats h96json:", h96["stats"], " stats val96 par:", {k:v for k,v in val["stats_par"].items() if k!='ninner'})
# and against h384 (should differ at h96 quadrature-error scale)
h384 = json.load(open("<project>/research/08_theory/K4e_tetrakis_d3/modelA_h384.json"))
print("val96 parallel vs modelA_h384.json (finer grid):")
for k in range(5):
    a = mp.mpf(val["Spar"][k]); b = mp.mpf(h384["Sdk"][k])
    print("  k=%d  |diff| = %s" % (k, mp.nstr(abs(a-b),3)))
