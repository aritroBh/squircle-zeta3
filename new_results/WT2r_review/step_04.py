base = '<project>/research/08_theory'
with open(base + '/WT2b_stabilization/result.md') as f:
    content = f.read()
print("LENGTH WT2b:", len(content))
with open(base + '/WT2_weight2_all_k/result.md') as f:
    c2 = f.read()
print("LENGTH WT2:", len(c2))
print(c2[:8000])
