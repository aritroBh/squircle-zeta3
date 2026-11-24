import json
p = '<scratch>'
data = json.load(open(p))
print(f"{'N1':>14} {'N2':>14} {'C':>14} {'ratio':>7}")
prev = None
for c in data:
    ratio = c[1]/c[0]
    print(f"{c[0]:>14.1f} {c[1]:>14.1f} {c[2]:>14.1f} {ratio:>7.3f}")
# check contiguity
ok = all(abs(data[i+1][0]-data[i][1]) < 1e-6*data[i][1] for i in range(len(data)-1))
print("contiguous:", ok, "| first N1:", data[0][0], "| last N2:", data[-1][1])
print("SUP C on [400,1e10]:", max(c[2] for c in data))
