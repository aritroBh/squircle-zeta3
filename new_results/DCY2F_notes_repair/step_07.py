base = "<project>/research/08_theory"
with open(f"{base}/DCY2R_referee_varpi1/result.md") as f:
    content = f.read()
# Print the middle part (lemma verdicts table) that was truncated
print(content[3000:6500])
