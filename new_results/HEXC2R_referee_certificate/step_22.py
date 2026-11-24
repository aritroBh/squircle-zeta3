
import os, time
os.makedirs("steps", exist_ok=True)
# reset progress log for this referee run
with open("progress.log", "w") as f:
    f.write(f"HEXC2R referee run started {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
print(open("../HEXC2_level6_certificate/progress.log").read()[-800:])
