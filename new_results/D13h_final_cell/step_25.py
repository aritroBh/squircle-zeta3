import json
d48 = json.load(open('final_cell.json'))
d192 = json.load(open('final_cell_J192.json'))
print("J48 :", {k: (v[:12] if isinstance(v,str) else v) for k,v in d48.items()})
print()
print("J192:", {k: (v[:12] if isinstance(v,str) else v) for k,v in d192.items()})
# sup from [100,400]
t = json.load(open('tight_100_400.json'))
sup1 = max(c[2] for c in t); w1 = [c[:2] for c in t if c[2]==sup1][0]
print("\nS1:", sup1, "at", w1)
s = json.load(open('sweep_cells_400_1e10.json'))
sup2 = max(c[2] for c in s); w2 = [c[:2] for c in s if c[2]==sup2][0]
print("S2:", sup2, "at", w2)
