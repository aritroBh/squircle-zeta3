txt = open('../DCY2_decay_first_correction/proof.md').read()
# print section 1 fully (from marker "## 1." to "## 2.")
i1 = txt.find("## 1."); i2 = txt.find("## 2.")
print(txt[i1:i2])
