import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
import glob
steps = sorted(glob.glob("step_*.py"))
print(steps[-6:])
print("="*80)
print(open(steps[-2]).read()[:2500])
