import subprocess, os
base = '<project>/research/08_theory/'
d3 = base+'TOW3_tower_depth6/'
# grep step files for which produce the pickles
out = subprocess.run(['grep','-l','L_parsed','/dev/null'] + [d3+f for f in os.listdir(d3) if f.startswith('step')], capture_output=True, text=True)
print('L_parsed in:', out.stdout)
out = subprocess.run(['grep','-l','e61_closed','/dev/null'] + [d3+f for f in os.listdir(d3) if f.startswith('step')], capture_output=True, text=True)
print('e61_closed in:', out.stdout)
out = subprocess.run(['grep','-l','arho56_simplified','/dev/null'] + [d3+f for f in os.listdir(d3) if f.startswith('step')], capture_output=True, text=True)
print('arho56_simplified in:', out.stdout)
