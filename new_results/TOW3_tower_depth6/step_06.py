base = '<project>/research/08_theory'
src = open(f'{base}/TOW2_tower_all_depths/tow2lib.py').read()
# print the truncated middle section
start = src.find('def complete_bell')
end = src.find('# ---------------- simplified tau')
print(src[start:end])
