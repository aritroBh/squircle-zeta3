content = open("<project>/research/08_theory/LEDGER.md").read()
idx = content.find("Corner angles")
print(content[idx-200:idx+4000])
