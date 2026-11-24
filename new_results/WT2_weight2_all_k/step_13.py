import os
for name in sorted(os.listdir('.')):
    print("="*30, name, "="*30)
    with open(name) as f:
        print(f.read()[:3000])
