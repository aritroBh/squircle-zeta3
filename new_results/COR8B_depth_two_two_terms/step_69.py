
import json, pickle
sy = json.load(open('../COR8_depth_two_any_angle/symbols.json'))
print(type(sy), list(sy.keys())[:20] if isinstance(sy, dict) else len(sy))
if isinstance(sy, dict):
    for k in list(sy.keys())[:5]:
        print("KEY:", k, "->", str(sy[k])[:300])
