import os, pickle
base = '<project>/research/08_theory/'
# Check engine dependencies
for p in ['OCT6_octagon_closed/oct5/Jloc.pkl', 'K1d_engine_fixed/Jtab_t.pkl']:
    fp = os.path.join(base,p)
    print(p, "exists:", os.path.exists(fp), os.path.getsize(fp) if os.path.exists(fp) else '')
# DOD1 assembly steps 34-38
for i in [34,35,36,37,38]:
    fn = base+f'DOD1_dodecagon/step_{i}.py'
    print("="*25, f"DOD1 step_{i}", "="*25)
    print(open(fn).read()[:1900])
