import json
base = "<project>/research/"
d = json.load(open(base + "08_theory/M4c_w5_closed_form/symbol.json"))
print("alphabet:", d["alphabet"])
print("delta32_symbol_pairing:", d["delta32_symbol_pairing"])
print("gpl_terms:", d["gpl_terms"], "controls:", d["product_and_Li5_controls"])
