import re
SRC = 'research/08_theory/K4c_3d_edge_law/'
tr = open(SRC+'transcript.json').read()
i = tr.find('h96outer (1194s)')
print(tr[i-50:i+1600].replace('\\n','\n').replace('\\"','"'))
