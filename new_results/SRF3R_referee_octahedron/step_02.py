import os
cwd = os.getcwd()
print("CWD:", cwd)
print(os.listdir(cwd))
print("--- parent ---")
print(os.listdir(os.path.dirname(cwd)))
print("--- grandparent ---")
print(os.listdir(os.path.dirname(os.path.dirname(cwd))))
