base = "<project>/research/08_theory"
with open(f"{base}/DCY2_decay_first_correction/proof.md") as f:
    content = f.read()
# find the §1 header
idx = content.find("## 1. Setup")
print(content[idx:idx+6500])
