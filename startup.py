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

def read_input(fn_params):
  # First read section 'Simulation_Parameters'
  with open(fn_params,"r") as f:
    f_list = [line.strip() for line in f.readlines()]
    for i,line in enumerate(f_list):
       if line.strip().startswith("SECTION SIMULATION_PARAMETERS"):
         continue_section = True
         section_counter = 0
         while continue_section:
           next_line = f_list[i+section_counter+1]
           if next_line.startswith("number_of_molecules"):
             glb.number_of_molecules = int(next_line.split()[-1])
           elif next_line.startswith("ncell"):
             glb.ncell = int(next_line.split()[-1])
           elif next_line.startswith("number_of_cycles"):
             glb.number_of_cycles = int(next_line.split()[-1])
           elif next_line.startswith("restart"):
             rs = next_line.split()[-1]
             if rs in ["TRUE","True","TRue","TRUe","true","t","T",".t.",".T."]:
               glb.restart = True
             else:
               glb.restart = False
           elif next_line.startswith("continue_run"):
             cr = next_line.split()[-1]
             if cr in ["TRUE","True","TRue","TRUe","true","t","T",".t.",".T."]:
               glb.continue_run = True
             else:
               glb.continue_run = False
           elif next_line.startswith("lpbc"):
             pb = next_line.split()[-1]
             if pb in ["TRUE","True","TRue","TRUe","true","t","T",".t.",".T."]:
               glb.lpbc = True
             else:
               glb.lpbc = False
           elif next_line.startswith("END SECTION SIMULATION_PARAMETERS"):
             continue_section = False
           else:
             print("Error reading Simulation_Parameters")

           # Should fail after 7, but either way, don't let infinite loop happen
           if section_counter > 100:
             print("Error: Check input file reader")

           section_counter += 1

  # Read section 'MOVES'
  with open(fn_params,"r") as f:
    f_list = [line.strip() for line in f.readlines()]
    for i,line in enumerate(f_list):
       if line.strip().startswith("SECTION MOVES"):

         continue_section = True
         section_counter = 0

         while continue_section:
           next_line = f_list[i+section_counter+1]
           if next_line.startswith("trans_prob"):
             glb.trans_prob = float(nlsplit[-1])
           elif next_line.startswith("trans_max_displ"):
             glb.trans_max_displ = float(nlsplit[-1])
           elif next_line.startswith("rot_max_displ"):
             glb.rot_max_displ = float(nlsplit[-1])
           elif next_line.startswith("drude_max_displ"):
             glb.drude_max_displ = float(nlsplit[-1])
           elif next_line.startswith("END SECTION MOVE_PROBS"):
             continue_section = False
           else: 
             print("Error reading Section Moves")
           if section_counter > 100:
             print("Check file reader")
             continue_section = False

           section_counter += 1

  # Read section 'SUPPLEMENTARY_INPUTS'
  with open(fn_params,"r") as f:
    f_list = [line.strip() for line in f.readlines()]
    for i,line in enumerate(f_list):
       if line.strip().startswith("SECTION SUPPLEMENTARY_INPUTS"):

         continue_section = True
         section_counter = 0

         while continue_section:
           next_line = f_list[i+section_counter+1]

           if next_line.startswith("box_file_names"):
             for b in range(glb.number_of_boxes):
               glb.box_file_names.append(next_line.split()[1+b])

           elif next_line.startswith("molecule_file_names"):
             for mt in range(glb.number_of_moltypes):
               glb.molecule_file_names.append(next_line.split()[1+mt])

           elif next_line.startswith("END SECTION SUPPLEMENTARY_INPUTS"):
             continue_section = False

           else:
             print("BAD JUST Bad codein go 2 BEd")

           section_counter += 1

           if section_counter > 100:
             print("How did you not hit an EOF, bro?")

  # Read section 'ELECTROSTATICS'
  with open(fn_params,"r") as f:
    f_list = [line.strip() for line in f.readlines()]
    for i,line in enumerate(f_list):
       if line.strip().startswith("SECTION ELECTROSTATICS"):

         continue_section = True
         section_counter = 0

         while continue_section:

           next_line = f_list[i+section_counter+1]

           if next_line.startswith("ewald_type"):
             glb.ewald_type = next_line.split()[-1]
           
           elif next_line.startswith("auto_ewald_par"):
             for b in range(glb.number_of_boxes):
               ae = next_line.split()[1+b]
               if ae in ["TRUE","True","TRue","TRUe","true","t","T",".t.",".T."]:
                 glb.auto_ewald_par.append(True)
               else:
                 glb.auto_ewald_par.append(False)

           elif next_line.startswith("ewald_conv_par"):
             for b in range(glb.number_of_boxes):
               cp = next_line.split()[1+b]
               if cp in ["TRUE","True","TRue","TRUe","true","t","T",".t.",".T."]:
                 glb.ewald_conv_par.append(True)
               else:
                 glb.ewald_conv_par.append(False)

           elif next_line.startswith("increase_cut_after"):
             for b in range(glb.number_of_boxes):
               ca = next_line.split()[1+b]
               if ca in ["TRUE","True","TRue","TRUe","true","t","T",".t.",".T."]:
                 glb.increase_cut_after.append(True)
               else:
                 glb.increase_cut_after.append(False)

           elif next_line.startswith("percent_cut_after"):
             for b in range(glb.number_of_boxes):
               glb.percent_cut_after.append(float(next_line.split()[1+b]))

           elif next_line.startswith("END SECTION ELECTROSTATICS"):
             continue_section = False

           else:
             print("Why didn't that cat come say hi to me?")
            
           section_counter += 1
   
           if section_counter > 100:
             print("How did you not hit an EOF, bro?")

def read_restart(fn_coords):
  with open(fn_coords,"r") as inp:
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

