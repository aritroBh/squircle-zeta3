import os, json
for f in sorted(os.listdir('.'), key=lambda s: (len(s), s)):
    p = os.path.join('.', f)
    print(f, os.path.getsize(p))
print()
# look at transcript: what was the last state?
try:
    t = json.load(open('transcript.json'))
    print(type(t), len(t) if hasattr(t, '__len__') else '')
except Exception as e:
    print("transcript load:", e)
