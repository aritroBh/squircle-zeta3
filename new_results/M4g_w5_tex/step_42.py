src=open('check_supplement.py').read()
i=src.find('# ---------------- (b)')
print(src[i:i+2600])
