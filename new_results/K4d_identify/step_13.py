import re
SRC = 'research/08_theory/K4c_3d_edge_law/'
tr = open(SRC+'transcript.json').read()
# find timing mentions like (123s)
for pat in ['h96outer', 'modelA', 'h64', 'h96 (', 'elapsed', '(', 's) stats']:
    pass
# Search for seconds patterns
mt = re.findall(r'\((\d+)s\)', tr)
print("elapsed seconds found:", mt[:40])
# search for 'h96' context with time
for m in re.finditer(r'h96', tr):
    i=m.start()
    seg = tr[i-60:i+80].replace('\\n',' ')
    if 's)' in seg or 'stats' in seg:
        print(seg)
