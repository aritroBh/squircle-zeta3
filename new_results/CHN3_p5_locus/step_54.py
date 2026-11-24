import time
time.sleep(240)
print(open('engdiag.log').read())
import os
print([f for f in ['eng12.pkl','vw12.pkl'] if os.path.exists(f)])
