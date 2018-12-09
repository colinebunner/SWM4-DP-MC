import numpy as np
import global_params as glb
import energy

def rigid_rot_move():

  rot_accept = 0

  # Store old energy for accepting/rejecting composite ANES rotation move
  energy_store = glb.box_energy
  energy_nok_store = glb.energy_nokspace

  # Pick random molecule for rotation
  trial_chain = int(np.floor(np.random.random()*glb.number_of_molecules))

  # Store original configuration
  xstore = np.copy(glb.xcoords)
  ystore = np.copy(glb.ycoords)
  zstore = np.copy(glb.zcoords)

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
  # when the move is accepted. Also, electronic moves are being performed with nested Markov
  # chain with modified potential in which the potential doesn't include k-space term.
  energy_elec_nok, loverlap = energy.sumup(lrecip=False)
  energy_elec_nok_store = energy_elec_nok
  energy_elec_full, loverlap = energy.sumup(lrecip=True)

  if loverlap:
    glb.xcoords = xstore
    glb.ycoords = ystore
    glb.zcoords = zstore

  else:
    # I don't think that the configurations will be too different if I ignore the k-space term in the Ewald sum, and so we will sample the
    # Drude positions with the scheme of Gelb [J. Chem. Phys. 118, 7747].
    nested_accept = False
    nestedx_store = np.copy(glb.xcoords)
    nestedy_store = np.copy(glb.ycoords)
    nestedz_store = np.copy(glb.zcoords)
    while not nested_accept:
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

        xd = glb.xcoords[trial_chain][4] - glb.xcoords[trial_chain][0]
        yd = glb.ycoords[trial_chain][4] - glb.ycoords[trial_chain][0]
        zd = glb.zcoords[trial_chain][4] - glb.zcoords[trial_chain][0]
        drudeDist = xd*xd + yd*yd + zd*zd

        if drudeDist > 0.04e0:
          glb.xcoords[trial_chain][4] -= dx_elec
          glb.ycoords[trial_chain][4] -= dy_elec
          glb.zcoords[trial_chain][4] -= dz_elec

        else:
          energy_new, loverlap = energy.sumup(lrecip=False)
    
          # Accept electronic move
          if np.random.random() < np.exp(-1.0*glb.beta_elec*(energy_new-energy_elec_nok)):
            energy_elec_nok = energy_new
          # Reject electronic move
          else:
            glb.xcoords[trial_chain][4] -= dx_elec
            glb.ycoords[trial_chain][4] -= dy_elec
            glb.zcoords[trial_chain][4] -= dz_elec

      # Accept/reject nested electronic move Markov chain
      energy_new_full, loverlap = energy.sumup(lrecip=True)
      deleprimed = energy_elec_nok-energy_elec_nok_store
      deletotal  = energy_new_full-energy_elec_full
      if np.random.random() < np.exp(-1.0*glb.beta_elec*(deletotal-deleprimed)):
        nested_accept = True
  
    # Accept composite move. Final energy_elec should be energy of system after all moves made.
    if np.random.random() < np.exp(-1.0e0*glb.beta*(energy_new_full-energy_store)):
      rot_accept += 1
      glb.box_energy = energy_new_full
      glb.energy_nokspace = energy_elec_nok
    else:
      glb.xcoords = xstore
      glb.ycoords = ystore
      glb.zcoords = zstore

  return rot_accept
