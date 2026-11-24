import os
base = '<project>'
# List the scripts directory
scripts = os.path.join(base, 'paper/squircle_note/scripts')
print(os.listdir(scripts))
print()
ieee = os.path.join(base, 'paper/ieee')
print(os.listdir(ieee))
