import os
base = os.path.dirname(os.getcwd())
with open(os.path.join(base,'OCT8_symbol_test','oct8.py')) as f:
    print(f.read())
