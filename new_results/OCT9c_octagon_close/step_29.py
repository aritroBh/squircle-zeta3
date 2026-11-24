import glob
for f in sorted(glob.glob('../OCT9_generic_octagon_closed/step_*.py')):
    src = open(f).read()
    if 'Ze_check' in src:
        print('====', f)
        i = src.find('Ze')
        print(src[:6000])
        break
