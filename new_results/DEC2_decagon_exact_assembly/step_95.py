base = "<project>/research/08_theory/OCT6_octagon_closed"
txt = open(f"{base}/step_61.py").read()
print(txt[:2500])
print("===== result.md excerpt =====")
md = open(f"{base}/result.md").read()
i = md.find("A-atom")
print(md[max(0,i-500):i+1800] if i!=-1 else "not found; search 'differentiat'")
