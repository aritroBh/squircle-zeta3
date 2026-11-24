src = open('../OCT9b_octagon_debug/edge_engine.py').read()
i = src.find('# ---------- atom constructors ----------')
j = src.find('# recursion cache')
print(src[i:i+4200])
