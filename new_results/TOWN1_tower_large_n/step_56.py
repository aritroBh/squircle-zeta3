import os
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
print("=== progress.log ===")
print(open('progress.log').read())
print("=== files ===")
for f in sorted(os.listdir('.')):
    if os.path.isfile(f):
        print(f"  {f} ({os.path.getsize(f)}b)")
print("  data/:", sorted(os.listdir('data')))
print("  results/:", sorted(os.listdir('results')))
