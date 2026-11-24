import pickle, subprocess
# quick nan check on the final itab140
IT = pickle.load(open("itab140.pkl","rb"))
nans = [m for m, v in IT.items() if 'nan' in v or '[' in v]
print("itab140: %d integrals, %d bad values" % (len(IT), len(nans)))
r = subprocess.run(["<project>/.venv/bin/python", "-u", "dod5d_assembly.py", "itab"],
                   capture_output=True, text=True, timeout=600)
print(r.stdout); print(r.stderr[-400:])
