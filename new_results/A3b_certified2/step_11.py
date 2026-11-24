
with open('../A3b_certified/result.md') as f:
    md = f.read()
print("total length:", len(md))
print(md[12000:28000])
