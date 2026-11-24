import os

base = "<project>/research/08_theory/"
brief = os.path.join(base, "COR3_classical_proof.brief")
with open(brief) as f:
    print(f.read())
