import re
SRC = 'research/08_theory/K4c_3d_edge_law/'
tr = open(SRC+'transcript.json').read()
# Find the h96outer block
i = tr.find('h96outer')
print("=== h96outer block ===")
print(tr[i-100:i+1400].replace('\\n','\n').replace('\\"','"'))
