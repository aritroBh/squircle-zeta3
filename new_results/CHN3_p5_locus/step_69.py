import time
time.sleep(300)
print(open('eng12b.log').read())
import os
print([f for f in ['eng12.pkl','vw12.pkl','forms1112.pkl'] if os.path.exists(f)])
