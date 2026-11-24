import os
def grep_files(roots, needles, ext='.md', maxcontent=0):
    hits=[]
    for root0 in roots:
        for root,dirs,files in os.walk(root0):
            for f in files:
                if f.endswith(ext):
                    p=os.path.join(root,f)
                    try: c=open(p,errors='ignore').read()
                    except: continue
                    if all(n in c for n in needles):
                        hits.append((p,len(c)))
    return hits
RES='research/08_theory/'
# Find result.md files discussing E2 / weight-2 with Li2 or dilog or "[d^2]"
for needles in [['Li_2','d^2'],['Li2','weight 2'],['dilog','edge'],['E2','Li']]:
    hits=grep_files([RES], needles)
    print(needles, "->", [os.path.basename(os.path.dirname(h[0]))+'/'+os.path.basename(h[0]) for h in hits[:6]])
