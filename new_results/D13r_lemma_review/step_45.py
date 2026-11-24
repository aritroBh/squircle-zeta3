import json
d='research/08_theory/D13h_final_cell'
cells = json.load(open(d+'/sweep_cells_400_1e10.json'))
print(type(cells), len(cells))
print(json.dumps(cells[0], indent=1)[:600])
print(json.dumps(cells[-1], indent=1)[:600])
