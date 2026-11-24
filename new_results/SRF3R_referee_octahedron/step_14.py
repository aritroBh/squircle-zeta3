import os
for root, dirs, files in os.walk("../../paper"):
    for f in files:
        if f.endswith(".tex"):
            print(os.path.join(root, f), os.path.getsize(os.path.join(root,f)))
