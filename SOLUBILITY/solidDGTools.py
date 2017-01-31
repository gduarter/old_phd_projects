#!/usr/bin/env python
import numpy as np
import sys

def calculate_density( nMol, gro_file ):
    # Open the original topology file and
    try:
        with open( gro_file, "r+") as gro:
            gro_lines = gro.readlines()
    except:
        raise IOError("No gro file was found")
    string = gro_lines[-1][3:-1]  # No boxes bigger than 1 nm are necessary here
    box = ' '.join(elem.split())
    box = box.split()
    if len(box) == 3:
        vol = float(box[0]) * float(box[1]) * float(box[2])
    else:
        v1, v2, v3 = np.array([float(box[0]),float(box[3]),float(box[4])]), np.array([float(box[1]),float(box[5]),float(box[6])]),np.array([float(box[2]),float(box[7]),float(box[8])])
        volume = np.dot(v1,v1) * np.dot(v2,v2) * np.dot(v3,v3) # volume in nm^3
    volume = volume * 1000. # volume in angstrom^3
    rho = nMol/volume 
    return rho

def kJ_to_kBT( energy_value, temperature ):
    """ Converts energy values from kJ/mol to kT

        VARIABLE          I/O         Type        Description
        ----------------------------------------------------------------------------
        energy_value      input       float       Energy in kJ/mol
        temperature       input       float       Temperature in kelvin
        new_energy_value  output      float       Energy in kBT

    """
    kB = 1.3806488 * 6.02214129/ 1000.0 # Boltzmann's constant (kJ/mol/K)
    beta = 1./(kB * temperature)
    new_energy_value = energy_value * beta
    return new_energy_value

def kBT_to_kJ( energy_value, temperature ):
    """ Converts energy values from kJ/mol to kT

        VARIABLE          I/O         Type        Description
        ----------------------------------------------------------------------------
        energy_value      input       float       Energy in kBT
        temperature       input       float       Temperature in kelvin
        new_energy_value  output      float       Energy in kJ/mol

    """
    kB = 1.3806488 * 6.02214129/ 1000.0 # Boltzmann's constant (kJ/mol/K)
    beta = 1./(kB * temperature)
    new_energy_value = energy_value / beta
    return new_energy_value

def convert_spring( kspring, temperature ):
    """ Calculates the spring constant to be used in Einstein Molecule simulations

        VARIABLE        I/O         Type        Description
        ----------------------------------------------------------------------------
        kspring         input       float       spring constant in kT/angstrom^2
        temperature     input       float       temperature
        kspring         output      float       spring constant in kJ/nm^2

    """
    return kBT_to_kJ(kspring * 100, temperature)

def modify_gro_file( gro_file, out_gro_file="solid.gro" ):
    """ renames the first molecule to FIX and the first atom to X1                                                         
        VARIABLE          I/O         Type        Description                                              
        ----------------------------------------------------------------------------                       
        gro_file          input       string      input gro file name                 
        out_gro_file      input       string      output gro file name                      
        nAtom             output      integer     number of atoms per molecule

    """
    try:
        with open( gro_file, "rb" ) as gro:
            gro_lines = gro.readlines()
    except:
        raise IOError("No coordinate file found")
    temp = []
    nAtom = 0
    for line in gro_lines:
        if ' 1UNK' in line:
            newline = line[0:3] + ' 1FIX' + line[8:]
            temp.append(newline)
            nAtom += 1
        else:
            temp.append(line)
    new_gro_lines = []
    for line in temp:
        if ' 1FIX' in line:
            if '1' in line.split():
                newline = line[0:13]+'X1'+line[15:]
                new_gro_lines.append(newline)
            else:
                new_gro_lines.append(line)
        else:
            new_gro_lines.append(line)
    with open( out_gro_file, "wb") as gro:
        gro.writelines( new_gro_lines )
    return nAtom

def thermal_wavelength( molar_mass ):
    """ Important: In order to avoid problems, I took in consideration that 
        kB = 1.38064852E-23 m^2 kg/s^2/ K and h = 6.626E-34 m^2 kg /s to simplify
        the calculation below. 
    """
    kB = 1.38064852 # m^2 kg/s^2/ K
    temperature = 298.15 # K
    PI = 3.14159265359
    hplanckSquared = 43.9048041749856E-45 #m^2 kg/ s
    avogadro = 6.02E23
    molar_mass= molar_mass/1000 # transform to kg/mol
    mass = molar_mass/avogadro

    return (np.sqrt(hplanckSquared/(2*PI*mass*kB*temperature)))*1.0E9 #convert to nanometers

def find_molar_mass(matrix, molecule):
    for elem in matrix:
        if molecule in elem:
            molar_mass = elem[3]
        else:
            print("Error: no molar mass")
            sys.exit()
    return molar_mass

def find_chemPot_solid_component( txtfile, nmols ):
    try:
        with open(txtfile, 'r') as table:
            for line in table.readlines():
                if 'TOTAL:' in line:
                    line = line.split()
                    DG = -float(line[-3]) # for MBAR
        return DG/nmols
    except IOError:
        return 'error'

def find_nmols( topfile ):
    with open( topfile ) as top:
        top_lines = top.readlines()
    for line in top_lines:
        if 'UNK' in line:
            nmols = float(line.split()[1]) + 1
    return nmols

def calculate_error_solid( txtfile, nmols ):
    with open( txtfile ) as txt:
        txt_lines = txt.readlines()
    count = []
    idx = 0
    for line in txt_lines:
        if 'States' in line:
            count.append[idx + 2]
        elif 'TOTAL:' in line:
            count.append[idx - 2]
        idx += 1
    errors = []
    for elem in txt_lines[count[0]:(count[1]+1)]:
        errors.append(float(elem.split()[-1]))
    errors = np.array(errors)/nmols
    return np.sqrt(np.dot(errors,errors))

def create_turnoff_forcefield( kspring, top_file, nAtom, ff_off_file ="morph.top" ):
    """ Organize 2 state topology for free energy calculation 
        where the force field is turned off.
        This script will work on files written by the tools we
        developed in the Mobley Lab.
        
        VARIABLE          I/O         Type        Description 
        ----------------------------------------------------------------------------  
        kspring           input       float       harmonic potential k_f in kT/(A^2)
        top_file          input       string      Initial topology file name
        nAtom             input       integer     number of atoms per molecule
        
    """
    # Open the original topology file and
    try:
        with open( top_file, "rb") as top:
            top_lines = top.readlines()
    except:
        raise IOError("No topology file was found")
    # extract key positions in the topology file
    idx = 0
    indices = []
    for line in top_lines:
        if '[ atomtypes ]' in line:
            indices.append([idx, '[ atomtypes ]'])
            idx += 1
        elif '[ moleculetype ]' in line:
            indices.append([idx,'[ moleculetype ]'])
            idx += 1
        elif '[ atoms ]' in line:
            indices.append([idx,'[ atoms ]'])
            idx += 1
        elif '[ bonds ]' in line:
            indices.append([idx,'[ bonds ]'])
            idx += 1
        elif '[ pairs ]' in line:
            indices.append([idx,'[ pairs ]'])
            idx += 1
        elif '[ angles ]' in line:
            indices.append([idx,'[ angles ]'])
            idx += 1
        elif '[ dihedrals ]' in line:
            indices.append([idx,'[ dihedrals ]'])
            idx += 1
        elif '[ system ]' in line:
            indices.append([idx,'[ system ]'])
            idx += 1
        elif '[ molecules ]' in line:
            indices.append([idx,'[ molecules ]'])
            idx += 1
        else: idx += 1
    # Modify atomtypes section and store it in list
    for elem in indices:
        if '[ atomtypes ]' in elem:
            count1 = elem[0]
            lim = count1 # To be used in the end
    for elem in indices:
        if '[ moleculetype ]' in elem:
            count2 = elem[0]
    atypes = top_lines[count1:count2-2]
    atypes.append('du             1    0.00000    1.008000  A             0             0\n')
    new_atomtypes = atypes
    # Create new moleculetype section
    for elem in indices:
        if '[ moleculetype ]' in elem:
            count1 = elem[0]
        if '[ atoms ]' in elem:
            count2 = elem[0]
    unk_moltype = top_lines[count1:count2-1]
    fix_moltype = top_lines[count1:count2-2]+['FIX          3\n']
    # Create new [ atoms ] section
    # All molecules
    for elem in indices:
        if '[ atoms ]' in elem:
            count1 = elem[0]
        if '[ bonds ]' in elem:
            count2 = elem[0]
    atoms = top_lines[count1+3:count2-1]
    final_unk = ['[ atoms ]\n']
    for elem in atoms:
        parameters = elem[:55] + '    12.0100      du   0.00000     12.0100\n'
        final_unk.append(parameters)
    # FIX molecule only
    final_fix = ['[ atoms ]\n']
    for elem in final_unk[1:]:
        if ' 1         ' in elem:
            newline = elem[:27]+'FIX     X1'+elem[37:]
            final_fix.append(newline)
        else:
            newline = elem[:27]+'FIX'+elem[30:]
            final_fix.append(newline)
    # Create [ position_restraint ] section
    pos_res_fix = ['[ position_restraints ]\n']
    i = 1
    while i < nAtom:
        i += 1                                                                                                                                                                                                               
        pos_res_fix.append('  %2d    1    %.2f %.2f %.2f %.2f %.2f %.2f\n' % (i, kspring, kspring, kspring, kspring, kspring, kspring))                                                                                      
    pos_res_unk = ['[ position_restraints ]\n']                                                                                                                                                                              
    j = 0                                                                                                                                                                                                                    
    while j < nAtom:                                                                                                                                                                                                         
        j += 1                                                                                                                                                                                                               
        pos_res_unk.append('  %2d    1    %.2f %.2f %.2f %.2f %.2f %.2f\n' % (j, kspring, kspring, kspring, kspring, kspring, kspring))                                                                                      
    # Build [ bonds ] section                                                                                                                                                                                                
    for elem in indices:                                                                                                                                                                                                     
        if '[ bonds ]' in elem:                                                                                                                                                                                              
            count1 = elem[0]                                                                                                                                                                                                 
        if '[ pairs ]' in elem:                                                                                                                                                                                              
            count2 = elem[0]                                                                                                                                                                                                 
    bond_list = top_lines[count1+3:count2-1]                                                                                                                                                                                  
    int_lines = []                                                                                                                                                                                                           
    for elem in bond_list:
        int_lines.append(elem.rstrip('\n').split() + [elem.rstrip('\n').split()[-2]] + [' 0.00000 \n']) #pay attention to the '\n' character
    final_bonds = ['[ bonds ]\n']
    for elem in int_lines:
        final_bonds.append(' '.join(elem))
    # Copy [ pairs ] section
    for elem in indices:
        if '[ pairs ]' in elem:
            count1 = elem[0]
        if '[ angles ]' in elem:
            count2 = elem[0]
    final_pairs = top_lines[count1:count2-1]
    # Build [ angles ] section
    for elem in indices:
        if '[ angles ]' in elem:
            count1 = elem[0]
        if '[ dihedrals ]' in elem:
            count2 = elem[0]
    angles_list = top_lines[count1+2:count2-1]
    int_lines = []
    for elem in angles_list:
        int_lines.append(elem.rstrip('\n').split() + [elem.rstrip('\n').split()[-2]] + [' 0.00000 \n'])
    final_angles = ['[ angles ]\n']
    for elem in int_lines:
        final_angles.append(' '.join(elem))
    # Build [ dihedrals ] section
    for elem in indices:
        if '[ dihedrals ]' in elem:
            count1 = elem[0]
        if '[ system ]' in elem:
            count2 = elem[0]
    dihedrals_list = top_lines[count1+2:count2-1]                                                                                                                                                                             
    int_lines = []                                                                                                                                                                                                           
    for elem in dihedrals_list:                                                                                                                                                                                              
        int_lines.append(elem.rstrip('\n').split() + [elem.rstrip('\n').split()[-3]] + [' 0.00000 '] + [elem.rstrip('\n').split()[-1]] + ['\n'])                                                                                     
    final_dihedrals = ['[ dihedrals ]\n']                                                                                                                                                                                                     
    for elem in int_lines:                                                                                                                                                                                                   
        final_dihedrals.append(' '.join(elem))                                                                                                                                                                               
    # Copy [ system ] section                                                                                                                                                                                                
    for elem in indices:                                                                                                                                                                                                     
        if '[ system ]' in elem:                                                                                                                                                                                             
            count1 = elem[0]                                                                                                                                                                                                 
        if '[ molecules ]' in elem:                                                                                                                                                                                          
            count2 = elem[0]                                                                                                                                                                                                 
    final_system = top_lines[count1:count2-1]                                                                                                                                                                                 
    # Build [ molecules ] section                                                                                                                                                                                            
    for elem in indices:                                                                                                                                                                                                     
        if '[ molecules ]' in elem:                                                                                                                                                                                          
            count = elem[0]                                                                                                                                                                                                      
    last = top_lines[count:]                                                                                                                                                                                                 
    nMol = int(last[-1][:-1].split()[1])                                                                                                                                                                                     
    fix, unk = 1, nMol - 1                                                                                                                                                                                                   
    fin = last[:2] + ['FIX       %s\n'%fix,'UNK       %s\n'%unk]                                                                                                                                                             
    # Build topology file list                                                                                                                                                                                               
    top1 = top_lines[:lim]+['\n']+new_atomtypes+['\n']+fix_moltype+['\n']+final_fix+['\n']+final_bonds+['\n']                                                                                                                 
    top2 = final_pairs+['\n']+final_angles+['\n']+final_dihedrals+['\n']+pos_res_fix +['\n']+unk_moltype+['\n']                                                                                                              
    top3 = final_unk+['\n']+final_bonds+['\n']+final_pairs+['\n']+final_angles+['\n']+final_dihedrals+['\n']+pos_res_unk+['\n']+fin                                                                                          
    new_top = top1+top2+top3                                                                                                                                                                                                 
    g = open(ff_off_file,"wb")
    for elem in new_top:
        g.write("%s" % elem)
    g.close()
    return

def create_turnoff_springs( kspring, nAtom, original_top, two_state_top="interacting.top", final_state_top="solid.top"):
    """ Gets information from both topologies and create a two state topology file
        VARIABLE          I/O         Type        Description 
        --------------------------------------------------------------------------------------------
        kspring           input       float       harmonic potential k_f in kT/(A^2)
        nAtom             input       integer     number of atoms per molecule
        original_top      input       string      name of the original topology for DA1 calculation
        two_state_top     input       string      name of the two-state topology for DA2 calculation
        final_state       input       string      name of the final state topology for DA1 calculation

    """
    try:
        with open(original_top) as top:
            top_lines = top.readlines()
    except:
        raise IOError("No topology file")
    idx = 0
    indices = []
    for line in top_lines:
        if '[ atomtypes ]' in line:
            indices.append([idx, '[ atomtypes ]'])
            idx += 1
        elif '[ moleculetype ]' in line:
            indices.append([idx,'[ moleculetype ]'])
            idx += 1
        elif '[ atoms ]' in line:                                                                                                                                                                                            
            indices.append([idx,'[ atoms ]'])                                                                                                                                                                                
            idx += 1                                                                                                                                                                                                         
        elif '[ bonds ]' in line:                                                                                                                                                                                            
            indices.append([idx,'[ bonds ]'])                                                                                                                                                                                
            idx += 1                                                                                                                                                                                                         
        elif '[ pairs ]' in line:                                                                                                                                                                                            
            indices.append([idx,'[ pairs ]'])                                                                                                                                                                                
            idx += 1                                                                                                                                                                                                         
        elif '[ angles ]' in line:                                                                                                                                                                                           
            indices.append([idx,'[ angles ]'])                                                                                                                                                                               
            idx += 1                                                                                                                                                                                                         
        elif '[ dihedrals ]' in line:                                                                                                                                                                                        
            indices.append([idx,'[ dihedrals ]'])                                                                                                                                                                            
            idx += 1                                                                                                                                                                                                         
        elif '[ system ]' in line:                                                                                                                                                                                           
            indices.append([idx,'[ system ]'])                                                                                                                                                                               
            idx += 1                                                                                                                                                                                                         
        elif '[ molecules ]' in line:                                                                                                                                                                                        
            indices.append([idx,'[ molecules ]'])                                                                                                                                                                            
            idx += 1
        else: idx += 1                                                                                                                                                                                                         
    # Create [ position_restraints ] section (two-state)                                                                                                                                                                     
    pos_res_fix = ['[ position_restraints ]\n']                                                                                                                                                                              
    i = 1                                                                                                                                                                                                                    
    while i < nAtom:                                                                                                                                                                                                         
        i += 1                                                                                                                                                                                                               
        pos_res_fix.append('  %2d    1    %.2f %.2f %.2f 0.00 0.00 0.00\n' % (i, kspring, kspring, kspring))                                                                                                                 
    pos_res_unk = ['[ position_restraints ]\n']                                                                                                                                                                              
    j = 0
    while j < nAtom:
        j += 1
        pos_res_unk.append('  %2d    1    %.2f %.2f %.2f 0.00 0.00 0.00\n' % (j, kspring, kspring, kspring))
    # Create FIX moleculetype section
    for elem in indices:
        if '[ moleculetype ]' in elem:
            count1 = elem[0]
        if '[ atoms ]' in elem:
            count2 = elem[0]
    fix_moltype = top_lines[count1:count2-2]+['FIX          3\n']
    # Create FIX [ atoms ] section 
    for elem in indices:
        if '[ atoms ]' in elem:
            count1 = elem[0]
        if '[ bonds ]' in elem:
            count2 = elem[0]
    atoms = top_lines[count1+3:count2-1]
    fix_atoms_new = ['[ atoms ]\n']
    for elem in atoms:
        if ' 1         ' in elem:
            newline = elem[:27]+'FIX     X1'+elem[37:]
            fix_atoms_new.append(newline)
        else:
            newline = elem[:27]+'FIX'+elem[30:]
            fix_atoms_new.append(newline)
    # Create [ molecules ] section
    for elem in indices:
        if '[ molecules ]' in elem:
            count = elem[0]
    last = top_lines[count:]                                                                                                                                                                                                 
    nMol = int(last[-1][:-1].split()[1])                                                                                                                                                                                     
    fix, unk = 1, nMol - 1                                                                                                                                                                                                   
    fin = last[:2] + ['FIX       %s\n'%fix,'UNK       %s\n'%unk]                                                                                                                                                             
    # Figure out place to put FIX topology and which parts of UNK will be repeated                                                                                                                                           
    for elem in indices:                                                                                                                                                                                                     
        if '[ moleculetype ]' in elem: Pos1 = elem[0]                                                                                                                                                                        
        if '[ bonds ]' in elem: Pos2 = elem[0]                                                                                                                                                                               
        if '[ molecules ]' in elem: Pos3 = elem[0]
        if '[ system ]' in elem: Pos4 = elem[0]                                                                                                                                                                           
    # Setting up UNK topology                                                                                                                                                                                                
    unk_top = top_lines[Pos1:Pos2] + pos_res_unk + ['\n'] + top_lines[Pos2:Pos3] + ['\n']                                                                                                                                    
    # Setting up FIX topology                                                                                                                                                                                                
    fix_top = fix_moltype +['\n']+ fix_atoms_new +['\n']+ pos_res_fix +['\n']+ top_lines[Pos2:Pos4] + ['\n']                                                                                                                 
    # Setting up two-state topology file                                                                                                                                                                                     
    new_top = top_lines[:Pos1] + fix_top + unk_top + fin                                                                                                                                                                     
    final_file = open(two_state_top,"wb")                                                                                                                                                                                    
    for elem in new_top:                                                                                                                                                                                                     
        final_file.write("%s" % elem)                                                                                                                                                                                        
    final_file.close()                                                                                                                                                                                                       
    return










