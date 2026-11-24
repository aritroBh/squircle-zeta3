base = "<project>/research/08_theory/"
import os
p = base+"LEDGER.md"
print(os.path.exists(p))
txt = open(p).read()
i = txt.rfind("DOD4")
print(txt[i-200:])
