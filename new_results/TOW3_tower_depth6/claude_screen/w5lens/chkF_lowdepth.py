import sys, os
sys.dont_write_bytecode=True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import prove_V56 as V
ok=True
for delta in (1,2,3,4):
    for s in range(0,delta//2+1):
        ok&=V.identity(delta,s)
print('low-depth non-adjacent identities with tow3lib forms, all True:',ok)
