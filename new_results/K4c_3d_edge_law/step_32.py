
import json
J = json.load(open("modelB_h32.json"))
for s in J["Sdk"]:
    print(repr(s[:80]))
