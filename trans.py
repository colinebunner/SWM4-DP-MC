import numpy as np
import global_params as glb
import energy

def trans_move():

  # Store old energy for accepting/rejecting composite ANES translation move
  energy_store = glb.box_energy

  # Pick random molecule for translation
  trial_chain = int(np.floor(np.random.random()*glb.number_of_molecules))

  # Store nuclear and electronic degrees of freedom
  xstore = glb.xcoords
  ystore = glb.ycoords
  zstore = glb.zcoords

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
  wrapped_x= False
  wrapped_y= False
  wrapped_z= False

  if xcom > glb.box_length:
    wrapped_x= True
    for i in range(5):
      glb.xcoords[trial_chain][i] = glb.xcoords[trial_chain][i]-glb.box_length
     
  if ycom > glb.box_length:
    wrapped_y= True
    for i in range(5): 
      glb.ycoords[trial_chain][i] = glb.ycoords[trial_chain][i]-glb.box_length

  if zcom > glb.box_length:
    wrapped_z= True
    for i in range(5): 
      glb.zcoords[trial_chain][i] = glb.zcoords[trial_chain][i]-glb.box_length

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

