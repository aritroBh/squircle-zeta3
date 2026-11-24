
import os
base = "<project>/"
print(os.listdir(base+"paper/squircle_note/scripts"))
print("="*60)
# find symeng.py
for root, dirs, files in os.walk(base+"research"):
    for f in files:
        if f == "symeng.py":
            print(os.path.join(root, f))
