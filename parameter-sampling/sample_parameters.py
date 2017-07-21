#!/usr/bin/env python
import numpy as np
import pylab as pl
import os, pickle, copy
import matplotlib.ticker as mtick
import random
from samplingtools import *
from smarty import *
from smarty import ForceField
import openeye
from openeye import oechem
from smarty.utils import get_data_filename
import argparse


def fake_loglikelihood():
    loglike = random.random()
    like = exp(loglike)
    return loglike, like


# Utility function to determine which parameters are in which molecules
def param_in_mols(labels, smirks):
    """Return list of True/False values as to whether specified SMIRKS pattern is in the molecules
    for which labels are provided. Labels should be as output by ForceField.labelMolecules"""
    smirks_in_mol = []
    for mol_entry in range(len(labels)):
        found = False
        for force in labels[mol_entry].keys():
            for (atom_indices, pid, s) in labels[mol_entry][force]:
                if s==smirks and not found:
                    smirks_in_mol.append(True)
                    found = True
        if not found:
            smirks_in_mol.append(False)
    return smirks_in_mol


def params_in_mols( labels, tuplelist ):
    """Return list of True/False values as to whether any of the specified SMIRKS patterns are in each of the specified molecules. Arguments are labels, as output by ForceField.labelMolecules, and smirkslist, list of SMIRKS of interest."""
    smirkslist = [tuplelist[i][1] for i in range(len(tuplelist))]
    params_is_in_mols = []
    i = 0
    for smirks in smirkslist:
        # Handle case where this is the first smirks
        if len(params_is_in_mols)==0:
            params_is_in_mols = param_in_mols(labels, smirks)
        # If this is a subsequent smirks, just update where needed
        else:
            tmp = param_in_mols(labels, smirks)
            # Make updates
            for idx in range(len(tmp)):
                if tmp[idx]:
                    params_is_in_mols[idx] = True
    return params_is_in_mols


# Define a uniform prior
def uniform_prior( theta, prange, parameter_key):
    """Defines a uniform prior

    Parameters
    ----------
    theta : float
        Parameter value to consider
    prange : dict
        Dictionary, keyed by parameter type, defining the parameter range to consider.
    parameter_key : str
        key from the dictionary which applies to this parameter
    """
    bounds = prange[parameter_key]
    if theta < bounds[0]:
        return 0
    elif theta > bounds[1]:
        return 0
    else: return 1


# Initialization function
def initialize(nmols, tuplelist, modify_key_by_type, ff, labels, prange, mean_unc, startparams):
    """Initialize for sampling.

    Parameters
    ----------
    nmols : int
        Number of molecules to use in sampling
    tuplelist : list of tuples
        SMIRKS strings of parameters we will sample
    modify_key_by_type : dict
        Dictionary, keyed by type, of parameters to modify
    ff : ForceField
        ForceField object to start with
    labels : list
        list of dictionaries with smirks labeled by type
    prange : dict
        Dictionary, keyed by parameter name, of what parameter values to consider
    mean_unc : float
        Uncertainty in reference energies
    startparams : str
        Choice of starting parameter values; use 'ff' to start with existing forcefield parameters; use 'uniform' to pick random starting point within range.

    Returns
    -------
    molidx_list : list of ints
        List of indexes of molecules we will use, of length nmols
    moveff : ForceField
        ForceField object we will be using for sampling
    measE : list of floats
        Measured energies for reference molecules
    meas_s_E : list of floats
        Uncertainties in measured energies for reference molecules
    start_param_vals : dict
        Dictionary of starting parameter values, keyed by SMIRKS and parameter key
    """

    # Find molecules for which this parameter changes uncertainty - list of True false falues
    params_is_in_mols = params_in_mols( labels, tuplelist )

    # Retrieve nmols molecules containing at least one of those params
    molidx_list = []
    ct=0
    while len(molidx_list) < nmols:# and ct < len(molidx_list):
        if params_is_in_mols[ct]:
            molidx_list.append(ct)
        ct+=1

    # Determine starting parameter values
    start_param_vals = {}
    for smirks_tuple in tuplelist:
        start_param_vals[smirks_tuple[1]] = {}
        for param_key in modify_key_by_type[smirks_tuple[0][0]]:
            thisrange = prange[param_key]
            if startparams=='uniform':
                start_param_vals[smirks_tuple[1]][param_key] = random.uniform( thisrange[0], thisrange[1])
            elif startparams=='ff':
                start_param_vals[smirks_tuple[1]][param_key] = ff.getParameter(smirks=smirks_tuple[1])[param_key]
            else:
                raise Exception("Error: Please choose starting parameter values.")

    # Retrieve measured and uncertainties
    measE = [ ref_energies[idx] for idx in molidx_list ]
    meas_s_E = np.array([mean_unc]*nmols)

    # Set starting parameter values and initialize force field
    moveff = copy.deepcopy(ff)

    # Initialize forcefield
    for smirks_tuple in tuplelist:
        params = moveff.getParameter(smirks=smirks_tuple[1])
        for param_key in modify_key_by_type[smirks_tuple[0][0]]:
            params[param_key]=str(start_param_vals[smirks_tuple[1]][param_key])
        moveff.setParameter(params, smirks=smirks_tuple[1])

    # Return stuff
    return molidx_list, moveff, measE, meas_s_E, start_param_vals


# Sampling function
def param_sampling(nsteps, nmols, oldE, tuplelist, modify_key_by_type, moveff,
                   move_range, prange, oemols, molidx_list, measE, meas_s_E,
                   start_param_vals):
    """ Monte Carlo sampling of parameters.
        
        Parameters
        ----------
        nsteps : int
        Number of MC steps
        nmols : int
        Number of molecules to use in sampling
        tuplelist : list of tuples
        tuples containing type and SMIRKS strings of parameters we will sample
        modify_key_by_type : dict
        Dictionary, keyed by type, of parameters to modify
        moveff : ForceField
        ForceField object we will be using for sampling
        prange : dict
        Dictionary, keyed by parameter name, of what parameter values to consider
        oemols : list of OEMols
        molidx_list : list of ints
        List of indexes of molecules we will use, of length nmols
        measE : list of floats
        Measured energies for reference molecules
        meas_s_E : list of floats
        Uncertainties in measured energies for reference molecules
        start_param_vals : dict
        Dictionary of starting parameter values, keyed by SMIRKS and parameter key
        
        Returns
        -------
        energies : numpy array
        Array with stored energy values
        param_vals : dict
        Dictionary of sampled parameters
        log_likelihood_t : numpy array
        Array containing the likelihoods of each step
        accepted_by_pkey : dict
        Dictionary containing accepted parameter values by key

        """

    # Init storage for energies and parameter values
    calcE = np.zeros((nmols), np.float32)
    energies = np.zeros((nsteps, nmols), np.float32)
    log_likelihood_t = np.zeros((nsteps), np.float32)
    param_vals = {}

    # create logfile
    outfile = open('%s_mol_%s_steps.log' % (nmols,nsteps),'w')

    # Evaluate old likelihood
    oldlog, oldlike = log_likelihood(oldE, measE, meas_s_E)
    accepted = 0
    #outfile.write("Old likelihood %.4g (log %.4f)\n" % (oldlike, oldlog))

    # Do MC of parameter moves
    accepted_by_pkey={}
    
    mcinfo = []
    
    for step in range(nsteps):
        stepinfo = []
        # Pick a random smirks to change
        indices = range(len(tuplelist))
        smirks_tuple = tuplelist[random.choice(indices)]
        # Create subdictionaries
        if smirks_tuple[0] not in accepted_by_pkey:
            accepted_by_pkey[smirks_tuple[0]] = {}
        # Pick a random parameter key to change
        param_key = random.choice(modify_key_by_type[smirks_tuple[0][0]])
        if not param_key in accepted_by_pkey[smirks_tuple[0]]:
            accepted_by_pkey[smirks_tuple[0]][param_key] = {}
            accepted_by_pkey[smirks_tuple[0]][param_key]['accepted']=0
            accepted_by_pkey[smirks_tuple[0]][param_key]['total']=0
        # Propose move
        chg = random.uniform(-move_range[param_key],move_range[param_key])
        # Get parameter value
        params = moveff.getParameter(smirks=smirks_tuple[1])
        paramval = float(params[param_key])
        newparam = paramval + chg

        # Focus only on things for which prior isn't zero
        while uniform_prior(newparam, prange, param_key)==0:
            chg = random.uniform(-move_range[param_key], move_range[param_key])
            newparam = paramval + chg
        
        # collect data for troubleshooting
        stepinfo.append(paramval)
        stepinfo.append(chg)
        stepinfo.append(newparam)
        
        # Prep for energy evaluation
        params[param_key]=str(newparam)
        moveff.setParameter(params, smirks=smirks_tuple[1])
        # Get new energy and likelihood
        for (idx, molid) in enumerate(molidx_list):
            # Pull molecule
            mol = oemols[molid]
            # Create system for energy eval
            topology = topologies[molid]
            system=moveff.createSystem(topology, [mol])
            # Get coordinates
            cpositions = reformat_oemol_coordinates(mol)
            # Get energy
            calcE[idx] = get_energy(system, cpositions[:,:,0])
            
            stepinfo.append(calcE[idx])
        
        newlog, newlike = log_likelihood(calcE, measE, meas_s_E)
        
        stepinfo.append(newlog)
        
        # Accept or reject move
        Pacc = np.exp(newlog - oldlog) * uniform_prior(newparam, prange, param_key)/uniform_prior(paramval, prange, param_key)
        accepted_by_pkey[smirks_tuple[0]][param_key]['total']+=1
        if step%500==0:
            # Print some progress info every N steps
            outfile.write("Step %s, Smirk: %s\n" % (step,smirks_tuple[1]))
            outfile.write("old %s: %.5f\n" % (param_key,paramval))
            outfile.write("change: %.5f\n" % chg)
            outfile.write("new %s: %.5f\n" % (param_key,newparam))
        accept = random.random() < Pacc
        if not accept:
            stepinfo.append("rejected")
            params[param_key] = str(paramval)
            moveff.setParameter(params, smirks=smirks_tuple[1])
        if accept:
            # Update "old" stuff if we accepted
            oldE = calcE
            oldlog, oldlike = newlog, newlike
            paramval = newparam
            accepted_by_pkey[smirks_tuple[0]][param_key]['accepted']+=1
            stepinfo.append("accepted")
        
        outfile.write("#steps accepted: %i\n" % accepted_by_pkey[smirks_tuple[0]][param_key]['accepted'])
        outfile.write("#steps sampled: %i\n" % accepted_by_pkey[smirks_tuple[0]][param_key]['total'])
        stepinfo.append(smirks_tuple[0])
        stepinfo.append(param_key)
        mcinfo.append(stepinfo)

        # Split parameters here as well.
        # Store data regardless
        energies[step][:] = oldE
        if step == 0:
            param_vals[smirks_tuple[1]] = {}
            param_vals[smirks_tuple[1]][param_key] = {}
        param_vals[smirks_tuple[1]][param_key][step] = paramval
        log_likelihood_t[step]=oldlog

        # Store values for any other parameters of this SMIRKS we're not changing
        for pkey in modify_key_by_type[smirks_tuple[0][0]]:
            if pkey != param_key:
                if step == 0:
                    param_vals[smirks_tuple[1]][pkey]={}
                    param_vals[smirks_tuple[1]][pkey][step] = float(params[pkey])
                if step > 0:
                    param_vals[smirks_tuple[1]][pkey][step] = param_vals[smirks_tuple[1]][pkey][step-1]

        # Store values for any parameters of other SMIRKS we're not changing
        for other_smirks_tuple in tuplelist:
            if other_smirks_tuple != smirks_tuple:
                params_tmp = moveff.getParameter(smirks=other_smirks_tuple[1])
                if step == 0:
                    param_vals[other_smirks_tuple[1]] = {}
                for pkey in modify_key_by_type[other_smirks_tuple[0][0]]:
                    if step == 0:
                        param_vals[other_smirks_tuple[1]][pkey] = {}
                        param_vals[other_smirks_tuple[1]][pkey][step] = float(params_tmp[pkey])
                    if step > 0:
                        param_vals[other_smirks_tuple[1]][pkey][step] = param_vals[other_smirks_tuple[1]][pkey][step-1]

    pickle.dump(mcinfo,open("MC_Info_%s_%s.p" % (nmols,nsteps),'wb'))

    for key in accepted_by_pkey.keys():
        for subkey in accepted_by_pkey[key].keys():
            print("Acceptance ratio of %s %s: %.2f" % (key, subkey, float(accepted_by_pkey[key][subkey]['accepted'])/float(accepted_by_pkey[key][subkey]['total'])))
            outfile.write("Acceptance ratio of %s %s: %.2f\n" % (key, subkey, float(accepted_by_pkey[key][subkey]['accepted'])/float(accepted_by_pkey[key][subkey]['total'])))

    outfile.close()
    

    return energies, param_vals, log_likelihood_t, accepted_by_pkey

## Main ##

# Input stuff
# Basic controls
parser = argparse.ArgumentParser(description='Samples forcefield parameters.')
parser.add_argument('-n','--nmols', type=int, help='number of molecules', required=True)
parser.add_argument('-ns','--nsteps', type=int, help='number of Monte Carlo steps', required=True)
parser.add_argument('-es','--equilsteps', type=int, help='number of equilibration steps',
                    required=False, default=0)
parser.add_argument('-t','--parmtype', type=str,
                    help='parameter type to be sampled; options: bond/angle/nonbonded',
                    required=False, default=None)

# Add number of molecules
nmols = parser.parse_args().nmols
nsteps = parser.parse_args().nsteps # MC steps to run
equilsteps = parser.parse_args().equilsteps # Discard data up to this step number from plots/stats
parmtype = parser.parse_args().parmtype # Parameter type to be sampled
if nmols == 0:
    print('Number of molecules needs to be greated than zero!\n')
    sys.exit()
elif nsteps <= equilsteps:
    print('Number of Monte Carlo steps needs to be bigger than the number of equilibration steps\n')
    sys.exit()

#nmols = 1
#nsteps = 1000
#equilsteps = 0
#parmtype = 'angle'

print("Number of molecules: %s" % nmols)
print("Number of Monte Carlo steps: %s" % nsteps)
print("Number of equilibration steps: %s" % equilsteps)
print("Parameter type: %s\n" % parmtype)


# Load pre-generated data if already generated, otherwise re-generate
if os.path.isfile('ref_energies.pickle') and os.path.isfile('topologies.pickle') and os.path.isfile('oemols_molnames.pickle'):
    file = open('ref_energies.pickle', 'r')
    ref_energies = pickle.load(file)
    file.close()
    file = open('topologies.pickle', 'r')
    topologies = pickle.load(file)
    file = open('oemols_molnames.pickle', 'r')
    oemols, mol_names = pickle.load(file)
    file.close()
    print("Loaded reference data.")
else:
    print("Re-generating reference data.")
    os.system('python generate_reference_data.py')
    file = open('ref_energies.pickle', 'r')
    ref_energies = pickle.load(file)
    file.close()
    file = open('topologies.pickle', 'r')
    topologies = pickle.load(file)
    file = open('oemols_molnames.pickle', 'r')
    oemols, mol_names = pickle.load(file)
    file.close()
    print("Loaded reference data.")

if os.path.isfile('s_E_mean.pickle'):
    file = open('s_E_mean.pickle', 'r')
    s_E, mean_unc = pickle.load(file)
    file.close()
else:
    print("Re-generating uncertainty data.")
    os.system('python get_uncertainties.py')
    file = open('s_E_mean.pickle', 'r')
    s_E, mean_unc = pickle.load(file)
    file.close()

# Take user inputs and apply functions
# What keys to modify for each smirks
modify_key_by_type = {'n':['rmin_half', 'epsilon'],'b':['length', 'k'], 'a':['angle','k'],'t':['k','periodicity','phase'], 'i':['k','periodicity','phase']}
# Figure out how to determine which molecules use the SMIRKS we're interested in
ffxml = get_data_filename('forcefield/Frosst_AlkEtOH.ffxml')
lb = ForceField(ffxml)
labels = lb.labelMolecules(oemols,verbose=False)
# Organize the types
ref_b = [label['HarmonicBondGenerator'] for label in labels]
#ref_t = [label['PeriodicTorsionGenerator'] for label in labels] #torsions are very complex
ref_n = [label['NonbondedGenerator'] for label in labels]
ref_a = [label['HarmonicAngleGenerator'] for label in labels]
#ref_labels = ref_b + ref_n + ref_a # + ref_t

if parmtype == 'bond':
    ref_labels = ref_b
elif parmtype == 'angle':
    ref_labels = ref_a
elif parmtype == 'nonbonded':
    ref_labels = ref_n
else:
    ref_labels = ref_b + ref_n + ref_a # + ref_t

# Extract list of reference SMIRKS types present in molecules, without repetition:
tuplelist = []
for param_set in ref_labels:
    if len(param_set) == 0: continue # monoatomic and torsionless molecules require this line
    else:
        for parameter in param_set:
            if (parameter[1],parameter[2]) not in tuplelist:
                tuplelist.append((parameter[1],parameter[2]))
            else: continue

# Do prep for parameter sampling
# Define parameter range to consider and move range
prange = {'rmin_half':[0.0, 3.], 'epsilon':[0.0, 1.0], 'length':[0.0,5.0], 'k':[0.,12000.], 'angle':[0.,180.]} #'k_{theta}':[10.,1000.]}#, 'k_{proper}':[0.,0.5], 'k_{improper}':[0.,20.], 'periodicity':[1.,5.], 'phase':[0.,180.]}
move_range = {'rmin_half':0.05, 'epsilon':0.01, 'length':0.005, 'k':1., 'angle':10.0}


# Initialize
molidx_list, moveff, measE, meas_s_E, start_param_vals = initialize(nmols, tuplelist, modify_key_by_type,
                                                                    lb, labels, prange, mean_unc, 'ff')

# Evaluate likelihood for our starting point
oldE = np.zeros( (nmols), np.float32 )
print("Evaluating initial likelihood...")
for (idx, molid) in enumerate(molidx_list):
    print("   Smiles string for this molecule %s" % OECreateIsoSmiString(oemols[idx]))
    # Pull molecule
    mol = oemols[molid]
    # Create system for energy eval
    topology = topologies[molid]
    system=moveff.createSystem(topology, [mol])
    # Get coordinates
    cpositions=reformat_oemol_coordinates(mol)
    # Get energy
    oldE[idx] = get_energy(system, cpositions[:,:,0])

#tuplelist = [('b0002','[#6X4:1]-[#1:2]')]
#tuplelist = [('a0002', '[#1:1]-[#6X4:2]-[#1:3]')]

print("Starting sampling procedure...")
# Do Monte Carlo Sampling
energies, param_vals, log_likelihood_t, accepted_by_pkey = param_sampling( nsteps, nmols, oldE, tuplelist,
                                                                          modify_key_by_type, moveff,
                                                                          move_range, prange, oemols,
                                                                          molidx_list, measE, meas_s_E,
                                                                          start_param_vals )

pickle.dump({'energies':energies, 'param_vals':param_vals, 'log_likelihood_t':log_likelihood_t, 'accepted_by_pkey':accepted_by_pkey}, open('results_%s_mol.p' % nmols,'wb'))


#nbins = 100
## Do stats/graphics for parameter sampling from above
#logfile = open("analysis_%s_%s_mol.log" % (parmtype,nmols), 'w')
#for element in tuplelist:
#    if parmtype == 'angle' and element[0][0] != 'a':
#        continue
#    if parmtype == 'bond' and element[0][0] != 'b':
#        continue
#    if parmtype == 'nonbonded' and element[0][0] != 'n':
#        continue
#    if parmtype == None: pass
#    # Plot histogram for each parameter under that smirks
#    print("Smirks %s:" % element[1])
#    for pkey in modify_key_by_type[element[0][0]]:
#        fig = pl.figure()
#        ax = fig.add_subplot(111)
#        dat = np.array([param_vals[element[1]][pkey][step] for step in range(len(param_vals[element[1]][pkey].keys()))][equilsteps:])
#        ax.hist(dat,bins=nbins,histtype='bar',edgecolor='black')
#        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.4f'))
#        ax.xaxis.set_major_locator(mtick.LinearLocator(numticks=5))
#        ax.yaxis.set_major_locator(mtick.LinearLocator(numticks=5))
#        ax.set_xlabel('Value for %s' % pkey)
#        ax.set_ylabel('Counts')
#        # Make title more compact
#        if len(element[1]) > 40:
#            title = element[1][:len(element[1])/2]+'\n'+element[1][len(element[1])/2:]
#        else:
#            title = element[1]
#        ax.set_title('%s' % title,fontsize=12)
#        ax.set_xlim((min(dat), max(dat)))
#        logfile.write("   Mean %s value %.5g (std %.5g)\n" % (pkey, dat.mean(), dat.std()))
#        logfile.write("   Force field parameter value %.5g\n" % float(lb.getParameter(element[1])[pkey]))
#        pl.tight_layout()
#        pl.savefig('histogram_%s_%s_%s.pdf' % (pkey, nmols, element[0]))
#        pl.close(fig)
#    # Plot parameter values on same plot with dual y axes
#    fig, ax1 = pl.subplots()
#    ax1.plot( np.arange(equilsteps,nsteps), [param_vals[element[1]][modify_key_by_type[element[0][0]][0]][step] for step in range(nsteps)][equilsteps:], 'b-')
#    ax1.set_xlabel('steps')
#    ax1.set_ylabel(modify_key_by_type[element[0][0]][0], color='b')
#    for tl in ax1.get_yticklabels():
#        tl.set_color('b')
#    ax2 = ax1.twinx()
#    ax2.plot( np.arange(equilsteps,nsteps), [param_vals[element[1]][modify_key_by_type[element[0][0]][1]][step] for step in range(nsteps)][equilsteps:], 'g-')
#    ax2.set_ylabel(modify_key_by_type[element[0][0]][1], color='g')
#    for tl in ax2.get_yticklabels():
#        tl.set_color('g')
#    # Make title more compact
#    if len(element[1]) > 40:
#        title = element[1][:len(element[1])/2]+'\n'+element[1][len(element[1])/2:]
#    else:
#        title = element[1]
#    ax.set_title('%s' % title,fontsize=12)
#    pl.tight_layout()
#    pl.savefig('pair_variation_%s_%s.pdf' % (nmols, element[0]))
#    pl.close(fig)
#    # Normalize data and reformat, then do plots
#    #probs_by_smirks={}
#    # Check if 2D parameter sampling; if so, grab data and do plots
#    if len(modify_key_by_type[element[0][0]])==2:
#        # Init storage
#        #probs_by_smirks[smirks] = np.zeros((nbins,nbins), np.float32)
#        # Pull data range
#        pkeys = modify_key_by_type[element[0][0]]
#        # Pull x and y values
#        xvals = [param_vals[element[1]][pkeys[0]][step] for step in range(nsteps)][equilsteps:]
#        yvals = [param_vals[element[1]][pkeys[1]][step] for step in range(nsteps)][equilsteps:]
#        # Plot
#        fig = pl.figure()
#        ax = fig.add_subplot(111)
#        #pl.hist2d(xvals, yvals, normed=True, norm=LogNorm())
#        hist = pl.hist2d(xvals, yvals, bins=40, normed=True, cmap='terrain')
#        pl.colorbar(hist[3],ax=ax)
#        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.4f'))
#        ax.xaxis.set_major_locator(mtick.LinearLocator(numticks=6))
#        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.4f'))
#        ax.yaxis.set_major_locator(mtick.LinearLocator(numticks=6))
#        ax.set_xlabel(pkeys[0])
#        ax.set_ylabel(pkeys[1])
#        # Make title more compact
#        if len(element[1]) > 30:
#            title = element[1][:len(element[1])/2]+'\n'+element[1][len(element[1])/2:]
#        else:
#            title = element[1]
#        ax.set_title('%s' % title,fontsize=12)
#        pl.tight_layout()
#        pl.savefig('probability_colorbar_%s_%s_%s_%s.pdf' % (nmols, element[0],
#                                                                  pkeys[0], pkeys[1]))
#        pl.close(fig)
#    # Plot likelihood vs step
#    fig = pl.figure()
#    ax = fig.add_subplot(111)
#    ax.plot(np.arange(equilsteps,nsteps), log_likelihood_t[equilsteps:])
#    # Make title more compact
#    if len(element[1]) > 30:
#        title = element[1][:len(element[1])/2]+'\n'+element[1][len(element[1])/2:]
#    else:
#        title = element[1]
#    ax.set_title('%s' % title,fontsize=12)
#    ax.set_xlabel('step')
#    ax.set_ylabel('log likelihood')
#    pl.savefig('log_likelihood_%s_%s.pdf' % (nmols,parmtype))
#    pl.close(fig)
#    # Print most likely parameter set found
#    mostlikelyval = log_likelihood_t.max()
#    idx = np.where(log_likelihood_t==mostlikelyval)
#    logfile.write("\nMost likely parameter values:\n")
#    for element in tuplelist:
#        line="   " + element[1] + ":  "
#        for st in idx[0]:
#            for pkey in modify_key_by_type[element[0][0]]:
#                line += "%s = %.3f, " % (pkey, param_vals[element[1]][pkey][st])
#            line+= "log likelihood %.4g" % log_likelihood_t[st]
#            logfile.write(line+"\n")
##    for st in idx[0]:
##        print("Energies at step %s:" % st, energies[st])
##    print("Reference energies:", measE)
#logfile.close()
#




