base='research/08_theory/A3c_claude_elementary/Q'
src=open(base+'/relations.py').read()
i=src.find('def small_calibration')
j=src.find('def main')
print(src[i:j])
