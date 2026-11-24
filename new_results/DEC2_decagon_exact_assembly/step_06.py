base = "<project>/research/08_theory"
with open(f"{base}/DEC1_decagon/result.md") as f:
    txt = f.read()
i = txt.find("D1 = ")
print(txt[i-200:i+2600])
