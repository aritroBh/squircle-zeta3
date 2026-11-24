
import json
fs = json.load(open('../COR4_M3alpha_classical/F_symbol.json'))
print(json.dumps(fs, indent=1)[:2500])
