print(open('brief.md').read()[:3000])
print('='*60)
for f in sorted(os.listdir('.')):
    if f.startswith('step_') or f=='transcript.json':
        print(f, os.path.getsize(f))
