
import os
os.chdir("<home>)")
print(os.getcwd())
# Read prev_result.md
try:
    with open("prev_result.md") as f:
        print(f.read())
except FileNotFoundError:
    print("prev_result.md NOT FOUND")
    # list files
    for fn in sorted(os.listdir(".")):
        print(fn)
