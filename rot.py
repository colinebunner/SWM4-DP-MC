import numpy as np
import global_params as glb
import energy

def rigid_rot_move():

  # Store old energy for accepting/rejecting composite ANES rotation move
  energy_store = glb.box_energy

  # Pick random molecule for rotation
  trial_chain = int(np.floor(np.random.random()*glb.number_of_molecules))

  # Store original configuration
  xstore = glb.xcoords
  ystore = glb.ycoords
  zstore = glb.zcoords

  # Rotate molecule around O atom
  alp = np.random.random()*2.0e0*np.pi
  bet = np.random.random()*2.0e0*np.pi
  gam = np.random.random()*2.0e0*np.pi
  ca = np.cos(alp)
  sa = np.sin(alp)
  cb = np.cos(bet)
  sb = np.sin(bet)
  cg = np.cos(gam)
  sg = np.sin(gam)

  xo = glb.xcoords[trial_chain][0]
  yo = glb.ycoords[trial_chain][0]
  zo = glb.zcoords[trial_chain][0]

  for iunit in range(1,5):
    # Move O atom to origin of coordinate system
    xu = glb.xcoords[trial_chain][iunit] - xo
    yu = glb.ycoords[trial_chain][iunit] - yo
    zu = glb.zcoords[trial_chain][iunit] - zo
    #Rotate and add position of O atom to coordinate
    glb.xcoords[trial_chain][iunit] = cg*(cb*xu+sa*sb*yu-ca*sb*zu) + sg*(ca*yu+sa*zu) + xo
    glb.ycoords[trial_chain][iunit] = sg*(-1.0e0*cb*xu-sa*sb*yu+sb*ca*zu) + cg*(ca*yu+sa*zu) + yo
    glb.zcoords[trial_chain][iunit] = sb*xu - sa*cb*yu + ca*cb*zu + zo

  # Don't want to recompute energy of old configuration every electronic move. Just update
  # when the move is accepted.
  energy_elec = energy.sumup()

  # Relec electronic moves. Pick a random drude oscillator to move.
  for em in range(glb.relec):
    trial_el = int(np.floor(np.random.random()*glb.number_of_molecules))

    # Electronic move will be simple cartesian displacement of Drude oscillators
    # Don't worry about wrapping because it's so unlikely that only the Drude oscillator
    # ends up outside of the box
    dx_elec = (np.random.random()-0.5e0)*glb.drude_max_displ
    dy_elec = (np.random.random()-0.5e0)*glb.drude_max_displ
    dz_elec = (np.random.random()-0.5e0)*glb.drude_max_displ

    glb.xcoords[trial_chain][4] += dx_elec
    glb.ycoords[trial_chain][4] += dy_elec
    glb.zcoords[trial_chain][4] += dz_elec

    energy_new = energy.sumup()

    # Accept electronic move
    if np.random.random() < np.exp(-1.0*glb.beta_elec*(energy_new-energy_elec)):
      energy_elec = energy_new
    # Reject electronic move
    else:
      glb.xcoords[trial_chain][4] -= dx_elec
      glb.ycoords[trial_chain][4] -= dy_elec
      glb.zcoords[trial_chain][4] -= dz_elec

  # Accept composite move. Final energy_elec should be energy of system after all moves made.
  if np.random.random() < np.exp(-1.0e0*glb.beta*(energy_elec-energy_store)):
    glb.box_energy = energy_elec
  else:
    glb.xcoords = xstore
    glb.ycoords = ystore
    glb.zcoords = zstore

