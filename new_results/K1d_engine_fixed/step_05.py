import os
for f in sorted(os.listdir('.')):
    sz = os.path.getsize(f)
    print(f"{sz:>10}  {f}")
