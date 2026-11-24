import os
pg5 = base+'PG5_alphabet_14gon/'
for f in sorted(os.listdir(pg5)):
    sz = os.path.getsize(pg5+f)
    print(f"{sz:>9}  {f}")
