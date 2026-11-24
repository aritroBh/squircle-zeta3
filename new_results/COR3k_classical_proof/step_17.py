import os
for f in ["step_01.py", "step_15.py"]:
    print(f"===== {f} =====")
    with open(f) as fh:
        print(fh.read()[:1500])
    print()
