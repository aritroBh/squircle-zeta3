
import os
print("CWD:", os.getcwd())
print("Contents:")
for d in sorted(os.listdir('.')):
    print(" ", d)
