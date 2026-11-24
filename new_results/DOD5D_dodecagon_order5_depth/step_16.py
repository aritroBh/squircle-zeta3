base = "<project>/research/08_theory/"
import glob
# The OCT5D step files are tiny auto-saved one-liners; the real code is in the transcript. Let's read a few key ones.
for n in [23,24,25,26]:
    fn = base+f"OCT5D_octagon_order5_depth/step_{n:02d}.py"
    print("="*100); print(fn)
    print(open(fn).read())
