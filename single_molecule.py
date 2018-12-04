import numpy as np
import global_params as glb
from startup import write_xyz
import energy
import matplotlib.pyplot as plt

# E-field strength (V/A)
efield = np.array([0.0e0,0.0e0,0.0e0])

steps = 50000
Relec = 10
istore_dipole = 20

# Allow rotation, or just check that oscillator stays at value of 0 [rot=False]
rot = False

# eA to Debye conversion (1 D = 0.20819434 eA)
dconv = 1.0e0/0.20819434e0

# eV to K conversion (1 eV = 11600 K) 
uconv = 11600.0e0

drude_max_displ = 0.04

telec     = 1.0e0
ttot      = 298.15e0
beta      = 1.0e0/ttot
beta_elec = 1.0e0/telec

#
def energy_extfield(mol):

  energy = 0.0e0

  # Drude portion
  drudex     = mol[4][0] - mol[0][0]
  drudey     = mol[4][1] - mol[0][1]
  drudez     = mol[4][2] - mol[0][2]
  drudeDist2 = drudex*drudex + drudey*drudey + drudez*drudez
  energy += 0.5e0*glb.kdrude*drudeDist2

  # External field portion (field in V/A, dipole is eA --> energy in eV: must conver to K)
  dipx, dipy, dipz = calc_dipole(mol)
  ux = dipx*efield[0]
  uy = dipy*efield[1]
  uz = dipz*efield[2]
  energy += (ux+uy+uz)*uconv
  
  return energy

def calc_dipole(mol):
  ### Induced dipole calculation ###
  # Don't forget: this doesn't correct for permanent dipole moment
  dipx = 0.0e0
  dipy = 0.0e0
  dipz = 0.0e0
  # Hopefully O bead counts for zero since we are rotating about O atom at origin
  for iunit in range(5):
    dipx += test_mol[iunit][0]*glb.qbeads[iunit]
    dipy += test_mol[iunit][1]*glb.qbeads[iunit]
    dipz += test_mol[iunit][2]*glb.qbeads[iunit]

  return dipx, dipy, dipz

def bin_prob_dist(dist,nbins):

  minval = np.min(dist)
  maxval = np.max(dist)
  #Round max and min values to nearest 0.01, add 0.01 so last value falls within a bin
  maxval = round(maxval,2)+0.01
  minval = round(minval,2)-0.01

  nobs = len(dist)

  bin_width = (maxval-minval)/nbins
  binned_dist = np.zeros(nbins)
  rbins = np.linspace(minval,maxval,num=nbins+1)
  rbins = np.array([(rbins[i+1]+rbins[i])/2. for i in range(nbins)])

  for el in dist:
    ibin = int((el-minval)//bin_width)
    binned_dist[ibin] += 1

  return rbins, binned_dist


if __name__ == "__main__":

  # Make SWM4-DP molecule, calculate initial energy 
  # Molecule is initially in x-y plane, dipole moment along
  # positive y axis
  test_mol = np.copy(glb.swm4dp_geom)
  energy_old = energy_extfield(test_mol)

  # For calculating average dipole moment
  dipole_tot = []
  dipole_x   = []
  dipole_y   = []
  dipole_z   = []

  for i in range(steps):
    if rot:

      # Copy is essential. I've been bit by this bug.
      geosave = np.copy(test_mol)

      # Rotate molecule around O atom
      alp = (np.random.random()-0.5e0)*glb.rot_max_displ
      bet = (np.random.random()-0.5e0)*glb.rot_max_displ
      gam = (np.random.random()-0.5e0)*glb.rot_max_displ
      ca = np.cos(alp)
      sa = np.sin(alp)
      cb = np.cos(bet)
      sb = np.sin(bet)
      cg = np.cos(gam)
      sg = np.sin(gam)

      xo = test_mol[0][0]
      yo = test_mol[0][1]
      zo = test_mol[0][2]

      for iunit in range(1,5):
        # Move O atom to origin of coordinate system
        xu = test_mol[iunit][0] - xo
        yu = test_mol[iunit][1] - yo
        zu = test_mol[iunit][2] - zo
        #Rotate and add position of O atom to coordinate
        glb.xcoords[trial_chain][iunit] = cg*(cb*xu+sa*sb*yu-ca*sb*zu) + sg*(ca*yu+sa*zu) + xo
        glb.ycoords[trial_chain][iunit] = sg*(-1.0e0*cb*xu-sa*sb*yu+sb*ca*zu) + cg*(ca*yu+sa*zu) + yo
        glb.zcoords[trial_chain][iunit] = sb*xu - sa*cb*yu + ca*cb*zu + zo
  
      energy_elec = energy_extfield(test_mol)
  
      # Relec electronic moves. Pick a random drude oscillator to move.
      for em in range(Relec):
        dx_elec = (np.random.random()-0.5e0)*drude_max_displ
        dy_elec = (np.random.random()-0.5e0)*drude_max_displ
        dz_elec = (np.random.random()-0.5e0)*drude_max_displ
  
        test_mol[4][0] += dx_elec
        test_mol[4][1] += dy_elec
        test_mol[4][2] += dz_elec
  
        energy_new = energy_extfield(test_mol) 
  
        # Accept electronic move
        if np.random.random() < np.exp(-1.0*beta_elec*(energy_new-energy_elec)):
          energy_elec = energy_new
        # Reject electronic move
        else:
          test_mol[4][0] -= dx_elec
          test_mol[4][1] -= dy_elec
          test_mol[4][2] -= dz_elec
  
      # Accept composite move. Final energy_elec should be energy of system after all moves made.
      if np.random.random() < np.exp(-1.0e0*beta*(energy_elec-energy_store)):
        energy_old = energy_elec
      else:
        test_mol = geosave 
  
    # Non-rotation case for calculating polarizability or just testing implementation
    else:

      geosave = np.copy(test_mol)
  
      energy_old = energy_extfield(geosave)
  
      dx_elec = (np.random.random()-0.5e0)*drude_max_displ
      dy_elec = (np.random.random()-0.5e0)*drude_max_displ
      dz_elec = (np.random.random()-0.5e0)*drude_max_displ
  
      test_mol[4][0] += dx_elec
      test_mol[4][1] += dy_elec
      test_mol[4][2] += dz_elec
  
      energy_new = energy_extfield(test_mol) 
  
      # Accept electronic move
      if np.random.random() < np.exp(-1.0*beta_elec*(energy_new-energy_old)):
        energy_old = energy_new
      # Reject electronic move
      else:
        test_mol[4][0] -= dx_elec
        test_mol[4][1] -= dy_elec
        test_mol[4][2] -= dz_elec

    if i%istore_dipole == 0:
      dipx, dipy, dipz = calc_dipole(test_mol)
      dipole_tot.append(np.sqrt(dipx*dipx+dipy*dipy+dipz*dipz))
      dipole_x.append(dipx)
      dipole_y.append(dipy)
      dipole_z.append(dipz)

  dipole_tot = np.array(dipole_tot)*dconv
  dtotBins, dtotDist = bin_prob_dist(dipole_tot,20)
  plt.plot(dtotBins,dtotDist)
  plt.show()
