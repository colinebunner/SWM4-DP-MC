import numpy as np
import matplotlib.pyplot as plt

# Maximum displacements of nuclear degrees of freedom
maxA = 1.0
maxB = 1.0

# Maximum displacements of electronic degrees of freedom
maxa = 0.1
maxb = 0.1

# How many electronic moves to attempt per nuclear move
Relec = 50

# How many simulation steps to perform
Nequil = 0
Nstep  = 200000
icollect = 10

nbins = 20

def nuc_move(A,B):

  # Pick a nuclear coordinate to move (equal probability for A and B)
  if np.random.random() < 0.5e0:
    A = A + (np.random.random()-0.5e0)*maxA
  else:
    B = B + (np.random.random()-0.5e0)*maxB

  return A, B

def elec_move(A,B,a,b,telec):
  enrgo = AaBb_energy(A,B,a,b)
  storea = a
  storeb = b
  # Pick an electronic coordinate to sample (equal prob  for a and b)
  if np.random.random() < 0.5e0:
    a = a + (np.random.random()-0.5e0)*maxa
  else:
    b = b + (np.random.random()-0.5e0)*maxb

  enrgn = AaBb_energy(A,B,a,b)

  # Accept
  if np.exp(((-1.0/telec)*(enrgn-enrgo))) > np.random.random():
    return a, b
  # Reject
  else:
    return storea, storeb

def AaBb_energy(A,B,a,b):
  return (0.5*(A**2.0e0 + B**2.0e0 - (A-B)**2.0e0 + (A-a)**2.0e0 + (A-b)**2.0e0 +
    (B-a)**2.0e0 + (B-b)**2.0e0 + (a-b)**2.0e0))


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

  return rbins

if __name__ == "__main__":

  posA =  0.0
  posB =  0.0
  posa =  0.0
  posb =  0.0

  distA = []
  distB = []
  dista = []
  distb = []

  # Simulation temperature and temperature of electronic d.o.f.
  temp_tot = 1.0
  temp_elec = 0.01


  for step in range(Nstep):

    energyOld = AaBb_energy(posA,posB,posa,posb)

    # Make a nuclear move
    trialA, trialB = nuc_move(posA,posB)
    # Call Relec electronic moves
    storea = posa
    storeb = posb
    for estep in range(Relec):
      posa, posb = elec_move(trialA,trialB,posa,posb,temp_elec)

    trialEnergy = AaBb_energy(trialA,trialB,posa,posb)

    # Final acceptance for total trial move
    if np.exp(((-1.0e0/temp_tot)*(trialEnergy-energyOld))) > np.random.random():
      posA = trialA
      posB = trialB
    else:
      posa = storea
      posb = storeb
#      print("Reject- Old: {:<8.5f}   New: {:<8.5f}".format(energyOld,trialEnergy))

    if step > Nequil and (step%icollect==0):
      distA.append(posA)
      distB.append(posB)
      dista.append(posa)
      distb.append(posb)
 
  # Bin into probability distribution
  Abins = bin_prob_dist(distA,nbins)
  Adist, Aedges = np.histogram(distA,bins=nbins,density=True)

  # Theoretical Gaussian probability distribution for A in Born-Oppenheimer limit
  Abins_theor = np.arange(min(Abins),max(Abins),0.01)
  Adist_theor = (1.0e0/(np.sqrt(temp_tot*2.0e0*np.pi)))*np.exp((-1.0e0*(Abins_theor**2.0e0))/(2.0e0*temp_tot))

  plt.plot(Abins,Adist,marker='o',mfc='white',linestyle='None')
  plt.plot(Abins_theor,Adist_theor,linestyle='-',color='black')
  plt.show()
  
