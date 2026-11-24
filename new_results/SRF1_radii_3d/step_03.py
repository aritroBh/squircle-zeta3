
# Set up my work directory and read the D16 engine + result.md
workdir = '<project>/research/08_theory/SRF1_radii'
os.makedirs(workdir, exist_ok=True)
print(open(d16 := '<project>/research/08_theory/D16_pure_weight/d16_engine.py').read())
