import numpy as np
import matplotlib.pyplot as plt

def oordf(xc,yc,zc,boxlen,gr,nrbins,rmin,rmax):
  rbw = float((rmax-rmin)/nrbins)
  xd = xc[1][0] - xc[0][0]
  yd = yc[1][0] - yc[0][0]
  zd = zc[1][0] - zc[0][0]
  xd = xd - boxlen*int(round(xd/boxlen))
  yd = yd - boxlen*int(round(yd/boxlen))
  zd = zd - boxlen*int(round(zd/boxlen))
  dist = np.sqrt(xd*xd+yd*yd+zd*zd)
  if dist < rmax:
    rbin = int((dist-rmin)//rbw)
    gr[rbin] += 2

def process_movie(fname):

  frame = 0
  nbins = 20
  rmin  = 2.5e0
  rmax  = 7.0e0
  groo  = np.zeros(nbins)
  groh  = np.zeros(nbins)
  rbw   = float((rmax-rmin)/nbins)
  rbins = np.array([((i+0.5)*rbw)+rmin for i in range(nbins)])
  boxl  = 14.0e0
  vol   = 14.0e0**3.0e0

  xcoords = np.zeros((2,5))
  ycoords = np.zeros((2,5))
  zcoords = np.zeros((2,5))

  with open(fname,"r") as f:
    for line in f:
      if line.startswith("10"):
        frame += 1
        beadcount = 0
      elif len(line.strip().split()) == 0:
        continue
      else:
        cl = line.strip().split()

        if beadcount < 5:
          molid = 0
          nbc = beadcount
        else:
          molid = 1
          nbc = beadcount -5
        
        xcoords[molid][nbc] = float(cl[1])
        ycoords[molid][nbc] = float(cl[2])
        zcoords[molid][nbc] = float(cl[3])
        beadcount += 1

      if beadcount == 9:
        oordf(xcoords,ycoords,zcoords,boxl,groo,nbins,rmin,rmax)
#        ohrdf(xcoords,ycoords,zcoords,boxl,groh,nbins,rmin,rmax)
  for rb in range(nbins):
    bvol = np.pi*(4.0e0/3.0e0)*(((rb+1)*rbw+rmin)**3.0e0-(rb*rbw+rmin)**3.0e0)*(2.0/vol)
    groo[rb] = (groo[rb]/bvol)*(1/2.0e0)
  return rbins, groo

rbins,oordf = process_movie("dimer-supert.xyz")
with open("groo.dat","w") as out:
  for i in range(len(rbins)):
    out.write("{:<8.3f} {:<8.3f}\n".format(rbins[i],oordf[i]))
