import os, json
p = '<scratch>'
print("exists:", os.path.exists(p))
if os.path.exists(p):
    data = json.load(open(p))
    print("num cells:", len(data))
    print("first:", data[:3])
    print("last:", data[-3:])
    print("max N2:", max(c[1] for c in data))
