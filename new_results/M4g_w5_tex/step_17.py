import json
w5=json.load(open('../M4e_w5_lchi/W5_5400.json'))
print(list(w5.keys()))
print(w5['value'][:120])
m4=open('../M4a_decomposition/M4_800.txt').read()
print(repr(m4[:200]))
