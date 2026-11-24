import os
base = '<project>'
for root, dirs, files in os.walk(os.path.join(base, 'research/08_theory')):
    level = root.replace(os.path.join(base,'research/08_theory'), '').count(os.sep)
    if level <= 2:
        indent = ' ' * 2 * level
        print('{}{}/'.format(indent, os.path.basename(root)))
        if level <= 1:
            subindent = ' ' * 2 * (level + 1)
            for f in files:
                print('{}{}'.format(subindent, f))
