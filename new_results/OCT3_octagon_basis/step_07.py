import json
from mpmath import mp, mpf, nstr
mp.dps = 280
d = json.load(open("consts62.json"))
print("name[0]:", d["names"][0])
print("val string:", d["vals"][0])
print("string length:", len(d["vals"][0]))
x0 = mpf(d["vals"][0])
print("parsed:", nstr(x0, 40))
# check the strings saved in step_04 more carefully
d31 = json.load(open("consts31.json"))
print("base[11] (Li3(1/2)):", d31["base"][11][:60], "... len", len(d31["base"][11]))
