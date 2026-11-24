import subprocess, glob
files = sorted(glob.glob('../OCT9_generic_octagon_closed/step_*.py'))
hits = [f for f in files if 'Ze' in open(f).read()]
print(hits)
src = open(hits[-1]).read()
i = src.find('Ze')
print(src[max(0,i-3000):i+2500])
