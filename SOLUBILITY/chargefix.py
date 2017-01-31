import parmed

mol2_filename = 'molecule.mol2'

#Correct the mol2 file partial atom charges to have a total net integer molecule charge                  
mol2f = parmed.formats.Mol2File
mol2f.write(parmed.load_file(mol2_filename).fix_charges(),mol2_filename)
