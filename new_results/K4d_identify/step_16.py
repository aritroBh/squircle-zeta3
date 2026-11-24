import re
SRC = 'research/08_theory/K4c_3d_edge_law/'
tr = open(SRC+'transcript.json').read()
i = tr.find('call_e7160b186707597a867f64445e4d7716', tr.find('h96outer'))
# print the tool result after this id
j = tr.find('call_e7160b186707597a867f64445e4d7716', i+10)
print(tr[j:j+1800].replace('\\n','\n').replace('\\"','"'))
