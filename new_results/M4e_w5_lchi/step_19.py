import json
p = '<project>/research/08_theory/M4d_w5_nice'
sym = json.load(open(p + '/nice_symbol.json'))
print("names:", sym['names'])
print("rank:", sym['rank'])
print("pivots:", sym['pivots'])
print("free:", sym['free'])
print("xp:", sym['xp'])
print("xq:", sym['xq'])
