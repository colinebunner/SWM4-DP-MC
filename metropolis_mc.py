import numpy as np
import global_params as glb
from startup import read_restart, init, write_xyz
from rot import rigid_rot_move
from trans import trans_move

def main():

  if glb.restart:
    read_restart(glb.restart_file)
  else:
    init()

  # Monte Carlo!
  for cycle in range(glb.number_of_cycles):
    for step in range(glb.number_of_molecules):
      # Pick what kind of move to perform
      move = np.random.random()
      if move < glb.trans_prob:
        trans_move()
      else:
        rigid_rot_move()
    write_xyz("test.xyz")

main()

