with open("proof.md") as f:
    txt = f.read()
old = """$$\nwith the sharper relative error $O(\\delta_N/\\Lambda)$ of [P3], step 4, valid as well; for
$\\delta=m$ (e.g. $f=\\sin\\frac{\\pi x}2$) this is $O(\\log n/n)$."""
# find the actual text
i = txt.find("with the sharper relative error")
print(repr(txt[i:i+220]))
