import os
print(open("progress.log").read())
print("FILES:", sorted(f for f in os.listdir('.') if not f.startswith('transcript')))
