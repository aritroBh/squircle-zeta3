import os
# Look at parent directories for referenced material
base = os.path.abspath('..')
print(base)
print(os.listdir(base))
