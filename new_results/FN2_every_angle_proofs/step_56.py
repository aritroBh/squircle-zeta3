with open("proof.md") as f:
    txt = f.read()
i = txt.find("T_3=\\frac{d_2}{a}")
print(txt[i-400:i+300])
