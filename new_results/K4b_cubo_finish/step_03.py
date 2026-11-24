d = "research/08_theory/K4b_cubo_finish"
import os
print(os.listdir(d))
with open(os.path.join(d, "prev_result.md")) as f:
    print(f.read())
