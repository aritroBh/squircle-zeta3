import sys; sys.dont_write_bytecode=True; sys.path.insert(0,'.')
src = open('hexc2_common.py').read()
i = src.find('def li_atom')
print(src[i:i+2200])
