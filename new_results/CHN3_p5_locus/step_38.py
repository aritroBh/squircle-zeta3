import time
time.sleep(60)
print(open('subres.log').read())
import os
print("files:", [f for f in ['subres1_A.txt','subres1_B.txt','subres0.txt'] if os.path.exists(f)])
