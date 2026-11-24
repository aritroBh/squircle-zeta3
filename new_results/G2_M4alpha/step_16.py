base = "<project>/research/"
import subprocess
# Look at verify_cor6.py and cor6_table.py to see how the pipeline was driven
with open(base + "08_theory/COR6_true_point_symbols/verify_cor6.py") as f:
    t = f.read()
print(t[:8000])
