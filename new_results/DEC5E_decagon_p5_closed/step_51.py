import os, time
for d in (1150, 1550):
    fn = f'build_{d}.log'
    if os.path.exists(fn):
        print('='*30, d)
        print(open(fn).read())
    else:
        print(d, "no log yet")
print([f for f in os.listdir('.') if f.startswith('gtab')])
