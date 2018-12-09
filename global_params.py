import numpy as np
### General simulation information ###
number_of_molecules = 64
ncell = 4
number_of_cycles = 50
restart = False
restart_file = "restart.xyz"
continue_run = False
# False for free dimer
lpbc = True

### Move type probabilities ###
trans_prob            = 0.500
rot_prob              = 1.000

### Max displacement stuff ###
trans_max_displ   = 0.4  # Angstrom
rot_max_displ     = 0.7  # Radians
drude_max_displ   = 0.04 # Angstrom

### Box info ###
# Cartesian boxlength in Angstrom (box is assumed cubic)
box_length = 12.417
# Temperature and reciprocal temperature in K and K^-1
temperature = 298.15e0
beta = 1.0e0/temperature
# Cutoff for nonbond interactions
rcut = 9.30 # Angstrom

### ANES-MC stuff ###
# Electronic temperature
temp_elec = 5.0e0 # K
beta_elec = np.float64(1.0e0/temp_elec)
# Number of electronic moves for ANES-MC
relec = 50

# Ewald sum stuff [definitely works now]
qtype = "Ewald"
nkvec = 7
kalp  = 3.2/rcut

# Output info
iwrite = 50
ienrg  = 20

# Bead information
nonbond_params = {101:[3.1589,93.2000,-1.77185],
 102:[0.0000,0.0000,0.55370],103:[0.0000,0.0000,-1.10740],
 104:[0.0e0,0.0e0,1.77185]}
bead_masses = {101:15.9990,102:1.0079,103:0.0000,104:0.0e0}
bead_names  = {101:"OW",102:"HW",103:"M",104:"D"}

#Fortran friendly form. 1st column is sigma, 2nd epsilon, 3rd charge, 4th beadtype
nbpar_fort = []
for i,bt in enumerate(nonbond_params.keys()):
  nbp = nonbond_params[bt]
  nbpar_fort.append([nbp[0],nbp[1],nbp[2],i+1])
nbpar_fort = np.array(nbpar_fort)

# Don't care I'm hardcoding this for calculating COM because 
masses = [15.9990,1.0079,1.0079,0.0e0,0.0e0,0.0e0]
# Ditto for Ewald sum
#qbeads = np.array([-1.77185,0.55370,0.55370,-1.10740,1.77185])
qbeads = np.array([0.0e0,0.55370,0.55370,-1.10740,0.0e0])

# SWM4-DP geometry
loh = 0.9572 # A
thetahalf = 52.26e0
angle = (90.0e0-thetahalf)*(np.pi/180.0e0)
xcomp = loh*np.cos(angle)
ycomp = loh*np.sin(angle)
swm4dp_geom = np.array([[0.0e0,0.0e0,0.0e0],
                        [xcomp,ycomp,0.0e0],
                        [-1.0e0*xcomp,ycomp,0.0e0],
                        [0.0e0,0.23808e0,0.0e0],
                        [0.0e0,0.0e0,0.0e0]])


### Information on each chain/unit ###
xcoords  = np.empty((number_of_molecules,5)) # 1st index is the molecule number, 
                                           # second is the unit. SWM4-DP has 5 units.
ycoords  = np.empty((number_of_molecules,5))
zcoords  = np.empty((number_of_molecules,5))

# Total energy of simulation box
box_energy = 0.0e0
energy_nokspace = 0.0e0

# Drude oscillator is tethered to the oxygen atom, which is bead 1, by a spring with the following force constant
kdrude = 503220.0 # K/A^2
qdrude = 1.77185  # units of e-1

# Tail correction to energy (only O has LJ params, so is only contributor to nonbonded tail correction)
sigrc  = 3.1589/rcut
sigrc3 = sigrc*sigrc*sigrc
sigrc9 = sigrc3*sigrc3*sigrc3
utailc = (8.0e0/3.0e0)*np.pi*(number_of_molecules/box_length**3.0e0)*nbpar_fort[0][1]*(nbpar_fort[0][0]**3.0e0)*((1.0e0/3.0e0)*sigrc9-sigrc3)*number_of_molecules
# False for finding dimer minimum energy
ltailc = True
