import numpy as np
import global_params as glb
from startup import read_restart, init, write_xyz
from rot import rigid_rot_move
from trans import trans_move
from energy import sumup

def main():

  if glb.restart:
    read_restart()
  else:
    init(rotate=True)

  # First energy calc
  glb.box_energy,loverlap = sumup()

  # Decide whether to continue run or overwrite
  write_xyz("run.xyz",(not glb.continue_run))
  # Monte Carlo!
  for cycle in range(glb.number_of_cycles):
    for step in range(glb.number_of_molecules):
      # Pick what kind of move to perform
      move = np.random.random()
      if move < glb.trans_prob:
        trans_move()
      else:
        rigid_rot_move()
    if cycle%glb.iwrite == 0:
      write_xyz("run.xyz",False)

  # Overwrite last restart file
  write_xyz("restart.xyz",True)

main()

