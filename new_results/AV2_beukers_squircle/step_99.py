import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
print(open('progress.log').read())
print("=== files ===")
print(sorted(os.listdir('.')))
import os.path
print("summary.md exists:", os.path.exists('summary.md'))
