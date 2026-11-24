# Check what words are in Lv
print("Lv words:", Lv.keys())

# Check which are missing from vals
for w in Lv.keys():
    if w not in vals:
        print(f"Missing: {w}")