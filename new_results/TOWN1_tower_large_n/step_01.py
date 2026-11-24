import os

base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')

# Create workdir
os.makedirs(workdir, exist_ok=True)
os.chdir(workdir)
print("Working dir:", os.getcwd())

# List the theory directories mentioned
for d in ['research/08_theory/WTS_weight_k_structure',
          'research/08_theory/TOW1_tower_law',
          'research/08_theory/WT2_weight2_all_k',
          'research/08_theory/WT3_weight3_all_k']:
    p = os.path.join(base, d)
    print("\n=== ", d)
    if os.path.exists(p):
        for root, dirs, files in os.walk(p):
            # limit depth
            depth = root.replace(p, '').count(os.sep)
            if depth <= 2:
                for f in files:
                    fp = os.path.join(root, f)
                    sz = os.path.getsize(fp)
                    print("  ", fp.replace(base+'/', ''), f"({sz}b)")
    else:
        print("  NOT FOUND")
