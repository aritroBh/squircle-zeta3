
import subprocess
# find thm:W4R in the paper
out = subprocess.run(["grep", "-n", "W4R", "../../../paper/ieee/main.tex"], capture_output=True, text=True)
print(out.stdout[:4000])
