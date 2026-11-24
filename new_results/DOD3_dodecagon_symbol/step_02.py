import os
parent = os.path.dirname(os.getcwd())
for d in os.listdir(parent):
    print(d)
