print(open("step_03.py").read()[:4000] if __import__('os').path.exists("step_03.py") else "none")
import json
d = json.load(open("transcript.json"))
print(type(d), len(d) if hasattr(d,'__len__') else '')
print(json.dumps(d, indent=1)[:1500])
