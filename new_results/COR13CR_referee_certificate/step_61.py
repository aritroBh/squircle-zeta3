import os
print(open('progress.log').read())
files = sorted(os.listdir('.'))
steps = sorted(f for f in files if f.startswith('step_'))
libs = [f for f in files if not f.startswith('step_')]
print("library/verifier files:", libs)
print("n step files:", len(steps))
