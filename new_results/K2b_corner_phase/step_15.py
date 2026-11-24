
import os
for dp,dn,fn in os.walk('research/08_theory'):
    for f in fn:
        if f=='result.md' and ('D13' in dp):
            print(os.path.join(dp,f))
