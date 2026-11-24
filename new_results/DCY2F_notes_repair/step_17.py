base = "<project>/research/08_theory"
with open(f"{base}/DCY2_decay_first_correction/result.md") as f:
    print(f.read())
print("#"*100)
with open(f"{base}/DCY2_decay_first_correction/result_nofinal.md") as f:
    print(f.read())
