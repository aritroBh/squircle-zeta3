import os, shutil, subprocess, time
wd = "<project>/research/08_theory/FN2_every_angle_proofs"
os.chdir(wd)
# remove stray dir created by my first exploratory block
if os.path.exists("research"):
    shutil.rmtree("research")
print(subprocess.run(["ls"],capture_output=True,text=True).stdout)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} cleaned stray dir; step files are auto-saved code blocks\n")
print("ok")
