import subprocess, glob, time
print("=== bgA768.log (current tail) ===")
print("\n".join(open("bgA768.log").read().splitlines()[-6:]))
print("ckpt files:", len(glob.glob("ckpt_h768/chunk_*.json")))
print("job pid:", open("jobA_h768.pid").read().strip(),
      "| alive:", subprocess.run(["pgrep","-f","jobA_h768.py"],capture_output=True,text=True).stdout.split())
print("worker children:", subprocess.run(["pgrep","-P",open("jobA_h768.pid").read().strip()],capture_output=True,text=True).stdout.split())
import os
print("\nfiles in work dir:", sorted(os.listdir(".")))
