
import os
print("CWD:", os.getcwd())
print("\nCWD contents:")
for f in sorted(os.listdir('.')):
    print("  ", f)
print("\n../D15b_polygon_law/ contents:")
for f in sorted(os.listdir('../D15b_polygon_law/')):
    print("  ", f)
