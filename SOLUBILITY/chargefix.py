#!/usr/bin/python
import parmed

mol2_filename = 'molecule.mol2'

mol2f = parmed.formats.Mol2File
mol2f.write(parmed.load_file(mol2_filename).fix_charges(),mol2_filename)


