import subprocess, os, time
r = subprocess.run(['uptime'], capture_output=True, text=True)
print(r.stdout)
# how long has eng12b been running, and is eng12.pkl appearing?
print(open('eng12b.log').read())
print([f for f in ['eng12.pkl','vw12.pkl','forms1112.pkl'] if os.path.exists(f)])
