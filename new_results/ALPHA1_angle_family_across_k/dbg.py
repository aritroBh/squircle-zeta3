
import sys
sys.path.insert(0, '.')
exec(open('step_30.py').read().split('Jc = to_hpl(Jexpr)')[0])
print("Jexpr terms:", len(Jexpr))
from collections import Counter
print("letters used:", Counter(a for (w,k) in Jexpr for a in w))
# try conv on each word, catch errors
fails = []
for (w, k), v in Jexpr.items():
    try:
        conv(w)
    except Exception as e:
        fails.append((w, str(e)))
print("conv failures:", len(fails))
for w, e in fails[:10]:
    print("   ", w, e)
# check a specific term: conv(('P','0','0')) etc
for w in [('P',), ('Q',), ('Q2',), ('P','0'), ('0','P'), ('1','Q2')]:
    try:
        print(w, "->", conv(w))
    except Exception as e:
        print(w, "ERR", e)
