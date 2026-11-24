import os
print(os.getcwd())
print(os.listdir("."))
print(os.listdir("..") if os.path.exists("..") else "no parent")
