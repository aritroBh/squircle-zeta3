import os
d = '../CHN3_p5_locus'
print("=== CHN3_p5_locus listing ===")
for f in sorted(os.listdir(d)):
    path = os.path.join(d, f)
    size = os.path.getsize(path) if os.path.isfile(path) else '<dir>'
    print(f"{size:>12}  {f}")
