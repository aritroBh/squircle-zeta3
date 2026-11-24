import os
# also look at LEDGER.md mentioned in the prompt
for root, dirs, files in os.walk("<project>"):
    for f in files:
        if "LEDGER" in f or "ledger" in f:
            print(os.path.join(root, f))
