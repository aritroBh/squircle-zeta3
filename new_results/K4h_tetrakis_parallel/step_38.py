import time
time.sleep(1500)
print(open("bgA768.log").read())
import subprocess
print("job alive:", subprocess.run(["pgrep","-f","jobA_h768.py"],capture_output=True,text=True).stdout.split())
