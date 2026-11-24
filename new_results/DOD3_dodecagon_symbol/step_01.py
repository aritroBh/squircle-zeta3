import os
base = os.getcwd()
print(base)
for d in os.listdir(base):
    print(d)
