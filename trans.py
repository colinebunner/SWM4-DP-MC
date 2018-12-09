import numpy as np
import global_params as glb
import energy

def trans_move():

  # Will return 0 if translation move is rejected, 1 if accepted
  trans_accept = 0

  # Store old energy for accepting/rejecting composite ANES translation move
  energy_store = glb.box_energy
  energy_nok_store = glb.energy_nokspace

  # Pick random molecule for translation
  trial_chain = int(np.floor(np.random.random()*glb.number_of_molecules))

  # Store nuclear and electronic degrees of freedom
  xstore = np.copy(glb.xcoords)
  ystore = np.copy(glb.ycoords)
  zstore = np.copy(glb.zcoords)

  # Random displacements for each cartesian direction
  dx_nuc = (np.random.random()-0.5e0)*glb.trans_max_displ
  dy_nuc = (np.random.random()-0.5e0)*glb.trans_max_displ
  dz_nuc = (np.random.random()-0.5e0)*glb.trans_max_displ
 
  # Add random displacement to all beads
  for i in range(5):
    glb.xcoords[trial_chain][i] += dx_nuc
    glb.ycoords[trial_chain][i] += dy_nuc
    glb.zcoords[trial_chain][i] += dz_nuc
  
  xcom = 0.0e0
  ycom = 0.0e0
  zcom = 0.0e0
  for i in range(5):
    # Get bead mass
    mass  = glb.masses[i]
    xcom  = xcom + glb.xcoords[trial_chain][i]*mass
    ycom  = ycom + glb.ycoords[trial_chain][i]*mass
    zcom  = zcom + glb.zcoords[trial_chain][i]*mass

  # More hardcoding. MW of water is 18.01528
  xcom = xcom/18.01528
  ycom = ycom/18.01528
  zcom = zcom/18.01528

  # Wrap if translation moves COM out of box in a direction. We won't let the
  # maximum translation get greater than half the box length, so we shouldn't
  # ever need to subtract more than boxlength from a coordinate.
  wrap_x_plus= False
  wrap_x_minus= False
  wrap_y_plus= False
  wrap_y_minus= False
  wrap_z_plus= False
  wrap_z_minus= False

  if xcom < 0.0e0:
    wrap_x_plus= True
  elif xcom > glb.box_length:
    wrap_x_minus= True

  if ycom < 0.0e0 :
    wrap_y_plus= True
  elif ycom > glb.box_length:
    wrap_y_minus= True

  if zcom < 0.0e0:
    wrap_z_plus= True
  elif zcom > glb.box_length:
    wrap_z_minus= True
  
  # Don't want to recompute energy of old configuration every electronic move. Just update
  # when the move is accepted. Also, electronic moves are being performed with nested Markov
  # chain with modified potential in which the potential doesn't include k-space term.
  energy_elec_nok, loverlap = energy.sumup(lrecip=False)
  energy_elec_nok_store = energy_elec_nok
  energy_elec_full, loverlap = energy.sumup(lrecip=True)

  # Kill move on hard-sphere overlap
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

        # Hard cutoff of 0.2 A on Drude displacements
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
      trans_accept += 1
      glb.box_energy = energy_new_full
      glb.energy_nokspace = energy_elec_nok

      # Wrap coordinates if trial move takes molecule outside of the box
      if wrap_x_plus:
        for iunit in range(5):
          glb.xcoords[trial_chain][iunit] += glb.box_length
      elif wrap_x_minus:
        for iunit in range(5):
          glb.xcoords[trial_chain][iunit] -= glb.box_length

      if wrap_y_plus:
        for iunit in range(5):
          glb.ycoords[trial_chain][iunit] += glb.box_length
      elif wrap_y_minus:
        for iunit in range(5):
          glb.ycoords[trial_chain][iunit] -= glb.box_length

      if wrap_z_plus:
        for iunit in range(5):
          glb.zcoords[trial_chain][iunit] += glb.box_length
      elif wrap_z_minus:
        for iunit in range(5):
          glb.zcoords[trial_chain][iunit] -= glb.box_length

    else:
      glb.xcoords = xstore
      glb.ycoords = ystore
      glb.zcoords = zstore


    return trans_accept
