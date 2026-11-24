base = "<project>/research/"
import subprocess
out = subprocess.run(["grep", "-rl", "F4_relation", base + "08_theory/COR5_corner_classicality/"], capture_output=True, text=True)
print(out.stdout)
out = subprocess.run(["grep", "-rl", "word_classes", base + "08_theory/COR5_corner_classicality/"], capture_output=True, text=True)
print(out.stdout)
