import time
time.sleep(720)
print(open("bgA768.log").read())
import subprocess, glob
print("ckpt files:", len(glob.glob("ckpt_h768/chunk_*.json")))
print("workers:", subprocess.run(["pgrep","-P","66142"],capture_output=True,text=True).stdout.split())
