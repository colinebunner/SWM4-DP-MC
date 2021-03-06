import numpy as np
import global_params as glb
import fenergy
from scipy.special import erf

kalp = glb.kalp
nkvec = glb.nkvec

# For converting sim units to K
kb   = np.float64(1.3806488E-23)
eps0 = np.float64(8.854187817E-12)
eu   = np.float64(1.602176565E-19)
qqfact = ((eu**2.0e0)*1.0e10)/(np.pi*4.0e0*eps0*kb)

# Self-interaction correction is constant in this simulation
qself = 0.0e0
for iunit in range(5):
  qself += glb.qbeads[iunit]**2.0e0
qself  = kalp*np.sqrt(1.0e0/np.pi)*glb.number_of_molecules*qself
qself *= qqfact

def sumup(lrecip=True):

  sys_energy = 0.0e0

  # Electrostatic energy (real + reciprocal - self - exclude)
  if glb.qtype == "Ewald":
    # Do this first because all beads have charge, so we can use electrostatic subroutine to get loverlap
    qreal,loverlap  = fenergy.energy.energy_qreal(glb.number_of_molecules,5,glb.box_length,glb.lpbc,glb.rcut,kalp,glb.xcoords,glb.ycoords,glb.zcoords,glb.qbeads)
    sys_energy += qreal
    sys_energy -= qself
    # Save time by not computing reciprocal space energy and off-diagonal exclusion terms
    if not loverlap:
      if lrecip:
        qrecip = fenergy.energy.energy_qkspace(glb.number_of_molecules,5,glb.box_length,kalp,nkvec,glb.xcoords,glb.ycoords,glb.zcoords,glb.qbeads)
        sys_energy += qrecip
      qexclude = fenergy.energy.energy_qexclude(glb.number_of_molecules,5,kalp,glb.xcoords,glb.ycoords,glb.zcoords,glb.qbeads)
      sys_energy -= qexclude
  else:
    qmimage,loverlap = fenergy.energy.energy_qmimage(glb.number_of_molecules,5,glb.box_length,glb.lpbc,glb.xcoords,glb.ycoords,glb.zcoords,glb.qbeads)
    sys_energy += qmimage

  # Kill move for overlap
  if loverlap:
    return sys_energy, True

  # Subroutine calculates Drude oscillator harmonic spring energy and nonbonded energy up to rcut
  enrg_nb = fenergy.energy.energy_nonbond(glb.number_of_molecules,5,glb.rcut,glb.box_length,glb.lpbc,glb.xcoords,
                     glb.ycoords,glb.zcoords,glb.nbpar_fort[0][0],glb.nbpar_fort[0][1])
  sys_energy += enrg_nb

  enrg_drude = fenergy.energy.energy_drude(glb.number_of_molecules,5,False,glb.xcoords,glb.ycoords,glb.zcoords,glb.kdrude)
  sys_energy += enrg_drude

  # Add tail correction (only O contributes to NB tail correction)
  if glb.ltailc:
    sys_energy += glb.utailc

  '''
  print("NB: {:<18.4f}".format(enrg_nb))
  print("NBTail: {:<18.4f}".format(glb.utailc))
  print("Drude: {:<18.4f}".format(enrg_drude))
  if glb.qtype == "Ewald" and lrecip:
    print("qreal: {:<18.4f}".format(qreal))
    print("qrecip: {:<18.4f}".format(qrecip))
    print("qself: {:<18.4f}".format(qself))
    print("qexclude: {:<18.4f}".format(qexclude))
    print("qs+qe: {:<18.4f}".format(qself+qexclude))
    print("qtot: {:<18.4f}".format(qreal+qrecip-qself-qexclude))
  elif not glb.qtype == "Ewald":
    print("qtot: {:<18.4f}".format(qmimage))
  print("---------------------------")
  print("Total: {:<18.4f}".format(sys_energy))
  '''

  return sys_energy, False
