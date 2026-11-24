import os
print("CWD:", os.getcwd())
print("CWD contents:")
for f in sorted(os.listdir('.')):
    print("  ", f, os.path.getsize(f) if os.path.isfile(f) else '<dir>')

# Test writability
try:
    with open('test_write.txt','w') as fh:
        fh.write('test')
    print("\nCan write to CWD: YES")
    os.remove('test_write.txt')
except Exception as e:
    print("\nCan write to CWD: NO -", e)
