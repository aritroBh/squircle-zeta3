with open("proof.md") as f:
    txt = f.read()
i = txt.find("## Theorem 2")
print(txt[i:i+7000])
