import subprocess
base = "<project>/research/08_theory/"
r = subprocess.run(["diff", base+"DOD4_dodecagon_closed/dod3_symbol.py", base+"DOD3_dodecagon_symbol/dod3_symbol.py"], capture_output=True, text=True)
print("dod3_symbol diff:", r.stdout[:2000] or "IDENTICAL")
r = subprocess.run(["diff", base+"DOD4_dodecagon_closed/dod3_R.py", base+"DOD3_dodecagon_symbol/dod3_R.py"], capture_output=True, text=True)
print("dod3_R diff:", r.stdout[:2000] or "IDENTICAL")
# look at key functions in dod3_symbol
txt = open(base+"DOD4_dodecagon_closed/dod3_symbol.py").read()
print(txt[:3500])
