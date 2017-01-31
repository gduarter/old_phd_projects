#!/usr/bin/env python
import solidDGTools as tools

old_kspring = 4000.
kspring = tools.convert_spring(old_kspring,temperature=298.15)
nAtom = tools.modify_gro_file('SYSTEM.gro')
tools.create_turnoff_forcefield( kspring, 'SYSTEM.top', nAtom )
tools.create_turnoff_springs( kspring, nAtom, 'SYSTEM.top' )


