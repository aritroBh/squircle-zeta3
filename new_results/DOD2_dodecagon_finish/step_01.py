import os
base = os.getcwd()
print("CWD:", base)
print("Listing CWD:")
for f in sorted(os.listdir('.')):
    print("  ", f)
print("\nListing parent:")
for f in sorted(os.listdir('..')):
    print("  ", f)
