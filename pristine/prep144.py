import ase.io
from ase.build import niggli_reduce
import numpy as np

unit_cell = ase.io.read('alpha_q.cif', index='0')
old_cell = unit_cell.get_cell()
unit_cell *= (2,1,1)
unit_cell.set_cell([
    old_cell[0] + old_cell[1],
    old_cell[0] - old_cell[1],
    old_cell[2]
])
print(unit_cell)

ase.io.write('test.xyz', unit_cell)
niggli_reduce(unit_cell)
print(unit_cell)
ase.io.write('sio2_ortho_cell.xyz', unit_cell)

defect_cell = unit_cell * (2,2,2)
ase.io.write("SiO_144.xyz",defect_cell)
del defect_cell[53]
ase.io.write('defect_144.xyz', defect_cell)
