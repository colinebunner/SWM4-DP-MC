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

  energy = []

  # Decide whether to continue run or overwrite
  write_xyz("dimer-min.xyz",(not glb.continue_run))
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
      write_xyz("dimer-min.xyz",False)
    if cycle%glb.ienrg == 0:
      etmp,lovr = sumup()
      energy.append(etmp)

  # Overwrite last restart file
  write_xyz("restart.xyz",True)

  if glb.number_of_cycles > 0:
    with open("dimer-min.out","w") as out:
      out.write("Simulation Information\n")
      out.write("Number of molecules: {}\n".format(glb.number_of_molecules))
      out.write("Temperature [K]: {}\n".format(glb.temperature))

      out.write("\nForce Field Information\n")
      out.write("Coulomb Calc. Type: {}\n".format(glb.qtype))
      if glb.qtype == "Ewald":
        out.write("num. k vectors: {}\n".format(glb.nkvec))
      out.write("Tail Corrections: {}\n".format(glb.ltailc))

      out.write("\nResults\n")
      out.write("Average Energy [K]: {:<16.4f} {:<16.4f}\n".format(np.mean(energy),np.std(energy,ddof=1)))

main()

