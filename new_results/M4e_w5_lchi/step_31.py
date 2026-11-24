import json
saved = json.load(open('<project>/research/08_theory/M4e_w5_lchi/W5_5400.json'))
print(list(saved.keys()))
print(len(saved['value']))
print(saved.get('bits'), saved.get('radius'))
print(saved['value'][:100])
