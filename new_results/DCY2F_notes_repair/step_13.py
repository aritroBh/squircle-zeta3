base = "<project>/research/08_theory"
with open(f"{base}/DCY2R_referee_varpi1/verify_dcy2r.py") as f:
    content = f.read()
i0 = content.find("# ---------------- 2. DISCREPANCY")
i1 = content.find("# ---------------- 5. region C")
print(content[i0:i1])
