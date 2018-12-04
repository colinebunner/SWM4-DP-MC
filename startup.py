import numpy as np
import global_params as glb

def init(rotate=False):

  atm = 0
  for i in range(glb.ncell):
    xi = (i+0.5) * (glb.box_length/glb.ncell)
    for j in range(glb.ncell):
      yj = (j+0.5) * (glb.box_length/glb.ncell) 
      for k in range(glb.ncell):
        zk = (k+0.5) * (glb.box_length/glb.ncell)
        if atm < glb.number_of_molecules:

          # Random displacements from regular lattice site to enhance equilibration
          rx = (np.random.random() - 0.5)*0.3
          ry = (np.random.random() - 0.5)*0.3
          rz = (np.random.random() - 0.5)*0.3
          if rotate:
            alp = np.random.random()*2.0e0*np.pi
            bet = np.random.random()*2.0e0*np.pi
            gam = np.random.random()*2.0e0*np.pi
            ca = np.cos(alp)
            sa = np.sin(alp)
            cb = np.cos(bet)
            sb = np.sin(bet)
            cg = np.cos(gam)
            sg = np.sin(gam)

          for iunit in range(5):
            x = glb.swm4dp_geom[iunit][0]
            y = glb.swm4dp_geom[iunit][1]
            z = glb.swm4dp_geom[iunit][2]
            if rotate:
              #Rotate 
              glb.xcoords[atm][iunit] = cg*(cb*x+sa*sb*y-ca*sb*z) + sg*(ca*y+sa*z) + xi + rx 
              glb.ycoords[atm][iunit] = sg*(-1.0e0*cb*x-sa*sb*y+sb*ca*z) + cg*(ca*y+sa*z) + yj + ry
              glb.zcoords[atm][iunit] = sb*x - sa*cb*y + ca*cb*z + zk + rz
            else:
              glb.xcoords[atm][iunit] = x + xi + rx 
              glb.ycoords[atm][iunit] = y + yj + ry
              glb.zcoords[atm][iunit] = z + zk + rz


          # Random displacement of Drude particle so Ewald sum doesn't break at beginning
          glb.xcoords[atm][4] += (np.random.random()-0.5e0)*0.01
          glb.ycoords[atm][4] += (np.random.random()-0.5e0)*0.01
          glb.zcoords[atm][4] += (np.random.random()-0.5e0)*0.01

          atm = atm + 1

def read_restart():
  with open(glb.restart_file,"r") as inp:
    f = [line.strip().replace("\n","") for line in inp.readlines()]
    for imol in range(glb.number_of_molecules):
      for iunit in range(5):
        cl = f[2+imol*5+iunit].split()
        glb.xcoords[imol][iunit] = float(cl[1])
        glb.ycoords[imol][iunit] = float(cl[2])
        glb.zcoords[imol][iunit] = float(cl[3])

def write_xyz(fn,overwrite=False):

  if overwrite:
    code = "w"
  else:
    code = "a+"

  with open(fn,code) as out:
    out.write("{:<8d}\n".format(5*glb.number_of_molecules))
    out.write("\n")
    for imol in range(glb.number_of_molecules):
      out.write("O {:<8.4f} {:<8.4f} {:<8.4f}\n".format(glb.xcoords[imol][0],glb.ycoords[imol][0],glb.zcoords[imol][0]))
      out.write("H {:<8.4f} {:<8.4f} {:<8.4f}\n".format(glb.xcoords[imol][1],glb.ycoords[imol][1],glb.zcoords[imol][1]))
      out.write("H {:<8.4f} {:<8.4f} {:<8.4f}\n".format(glb.xcoords[imol][2],glb.ycoords[imol][2],glb.zcoords[imol][2]))
      out.write("M {:<8.4f} {:<8.4f} {:<8.4f}\n".format(glb.xcoords[imol][3],glb.ycoords[imol][3],glb.zcoords[imol][3]))
      out.write("D {:<8.4f} {:<8.4f} {:<8.4f}\n".format(glb.xcoords[imol][4],glb.ycoords[imol][4],glb.zcoords[imol][4]))

