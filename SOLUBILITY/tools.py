# =============================================================================
# Adapted from David Mobley's script for PharmSci 272
#   and scripts process_submissions.py and uncertainty_check.py from SAMPL4
# By Caitlin C. Bannan
# February 2016
# Additions by Guilherme D. R. Matos
# January 2017
# =============================================================================

import numpy
from numpy import *
import scipy
import scipy.stats
import matplotlib. pyplot as plt
import matplotlib.patches as patches
from pylab import *
import scipy.integrate
from operator import itemgetter, attrgetter

def colors():
    '''Function that generates a random hexadecimal color code
        '''
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())

def rmserr( set1, set2 ):
    """Compute and return RMS error for two sets of equal length involving the same set of samples."""
    tot = 0.
    for i in range(len(set1)):
        tot+= (set1[i] - set2[i])**2
    return sqrt(tot/float(len(set1)))

def correl(x,y ):
    """For two data sets x and y of equal length, calculate and return r, the product moment correlation. """
    xar = numpy.array(x)
    yar = numpy.array(y)


    #Compute averages (try not to require numerical library)
    avgx=sum(xar)/float(len(xar))
    avgy=sum(yar)/float(len(yar))

    #Compute standard deviations
    sigmax_sq=0
    for elem in xar:
        sigmax_sq+=(elem-avgx)**2
    sigmax_sq=sigmax_sq/float(len(xar))
    sigmay_sq=0
    for elem in yar:
        sigmay_sq+=(elem-avgy)**2
    sigmay_sq=sigmay_sq/float(len(yar))

    sigmax=sqrt(sigmax_sq)
    sigmay=sqrt(sigmay_sq)

    #Compute numerator of r
    num=0
    for i in range(len(xar)):
        num+=(xar[i]-avgx)*(yar[i]-avgy)
    #Compute denominator of r
    denom=len(xar)*sigmax*sigmay

    corr = num/denom
    return corr

def percent_within_half( x, y, compare = 0.5):
    """
    input: x and y should be equal length arrays
    compare is a float the difference will be compared to
    returns: percent of values within the 'compare' value
    """

    diff = x - y
    indices = where( abs(diff) < compare)
    return 100.*float(len(indices[0]) )/float(len(x))    

def RemoveNones(data, remove=None):
    """
    Given a numpy array of arrays this will remove columns in the array where the value for any row is 'remove'
    input: data = numpy array  of any length with any number of rows
        remove = the value that if found in any row that entire column will be removed
    """
    shape = len(data.shape)
    if shape == 1:
        test  = array([data]) 
    else:
        test = data.copy()
    
    delRows = []
    for i in range(len(test[0])): 
        if None in list(test[:,i]):
            delRows.append(i)

    test = delete(test, delRows, 1)
    if shape == 1:
        return test[0]
    else:
        return test

def CorrectSign(x, y):
    """
    input: two arrays (or lists) of the same length
    output: percent of the values that agree with sign
    """
    xar = numpy.array(x)
    yar = numpy.array(y)
    agree = float(numpy.sum(numpy.sign(xar) == numpy.sign(yar)))
    return agree/float(len(x)) * 100.0

def MeanSignDeviation(x,y):
    """
    Calculates the mean signed deviation (mean sign error)
    that is MSD = sum(y_i - x_i) / n
    """
    x = array(x)
    y = array(y)
    sumDif = sum(y-x)
    
    return sumDif / len(x)
    
# Copied from David L. Mobley's scripts written for SAMPL4 analysis (added calculation uncertainty)
def bootstrap( datasets ):
    """Take a list of datasets of equal length. Use bootstrap to construct a new list of datasets of the same length and return them. The procedure is to select random data point indices (running from 0 to the set length, with replacement) and construct new datasets from those indices from each of the specified sets."""

    #Initialize storage; determine length of set
    newsets = []
    npoints = len( datasets[0] )
    
    #Pick random datapoint indices
    idx = numpy.random.randint( 0, npoints, npoints) #Create an array consisting of npoints indices, where each index runs from 0 up to npoints. 

    for dataset in datasets:
        newsets.append( dataset[idx] )
        #Error check
        if len(dataset) <> npoints:
            raise BaseException("Error: Variable length datasets passed to bootstrap function, which is not acceptable. Terminating.")

    return newsets 

# Copied from David L. Mobley's scripts written for SAMPL4 analysis (added calculation uncertainty)
def bootstrap_exptnoise( calc1, expt1, exptunc1, returnunc = False):
    """Take two datasets (equal length) of calculated and experimental values. Construct new datasets of equal length by picking, with replacement, a set of indices to use from both sets. Return the two new datasets. To take into account experimental uncertainties, random noise is added to the experimental set, distributed according to gaussians with variance taken from the experimental uncertainties. Approach suggested by J. Chodera.
Optionally, 'returnunc = True', which returns a third value -- experimental uncertainties corresponding to the data points actually used."""
    
    # Make everything an array just in case
    calc = array(calc1)
    expt = array(expt1)
    exptunc = array(exptunc1)
    npoints = len(calc)

    #Pick random datapoint indices
    idx = numpy.random.randint( 0, npoints, npoints) #Create an array consisting of npoints indices, where each index runs from 0 up to npoints.

    #Construct initial new datasets
    newcalc = calc[idx]
    newexpt = expt[idx]
    newuncExp = exptunc[idx]

    #Add noise to experimental set
    noise = numpy.random.normal( 0., exptunc) #For each data point, draw a random number from a normal distribution centered at 0, with standard devaitions given by exptunc
    newexpt += noise

    if not returnunc:
        return newcalc, newexpt
    else:
        return newcalc, newexpt, newuncExp

# Copied from David L. Mobley's scripts written for SAMPL4 analysis (added calculation uncertainty)
# CCB: I want to add an optional bootstrap with or without noise, I think this can be done without defining a separate function
def stats_array( calc1, expt1, exptunc1, boot_its, sid = "Number?", noise = True):
    """Compute statistics given sets of calculated and experimental values, experimental uncertainties, and a number of bootstrap iterations, for error estimates.
    Returns: 
        avgerr/d_avgerr: Average error and error estimate
        RMS/d_RMS: RMS error and error estimate
        AUE/d_AUE: AUE and error estimate
        tau/d_tau: Kendall tau and error estimate REMOVED
        R, d_R: Pearson R and error estimate
        maxE, d_maxE: Maximum error and uncertainty
        per, d_per: percent calculated with correct sign
        """
    # Assign array type to gaurentee calculations will work correctly
    # Allows users to submit anything "list like" 
    calc = np.array(calc1, dtype = np.float64)
    expt = np.array(expt1, dtype = np.float64)
    exptunc = np.array(exptunc1, dtype = np.float64)

    #COMPUTE OVERALL PERFORMANCE
    #Compute various metrics -- average error, average unsigned error, RMS error, Kendall Tau, Pearson correlation coefficient and Percent Correct Sign
    per = CorrectSign(calc, expt) 
    #Avg err
    avgerr = (calc - expt).mean()
    #RMS err
    RMS = numpy.sqrt( 1./float(len(expt)) * ((expt-calc)**2).sum() )
    #Average unsigned error
    AUE = numpy.abs( calc - expt).mean()
    #Tau
    tau, ptau = scipy.stats.kendalltau( expt, calc)
    #Pearson R
    R, pr = scipy.stats.pearsonr( expt, calc )
    #Max
    maxE = max(abs(calc-expt))
    r2 = correl(expt, calc) 
    #DO BOOTSTRAP ERROR ANALYSIS FOR OVERALL PERFORMANCE
    print "   Bootstrap for %s..." % sid
    avgerrs = numpy.zeros( boot_its)
    RMSerrs = numpy.zeros( boot_its)
    AUEs = numpy.zeros(boot_its)
    taus = numpy.zeros(boot_its)
    Rs = numpy.zeros(boot_its)
    maxEs = numpy.zeros(boot_its)
    pers = numpy.zeros(boot_its)
    r2s = numpy.zeros(boot_its)
    for it in range(boot_its):
        if noise:
            [ c, e ] = bootstrap_exptnoise( calc, expt, exptunc)  #Generate new sets, adding appropriate random variation to the experimental set.
        else:
            [ c, e ] = bootstrap( [calc, expt] ) # Generate new sets without adding noise to the experimental set
        avgerrs[it] = (c-e).mean()
        RMSerrs[it] =  numpy.sqrt( 1./float(len(e)) * ((e-c)**2).sum() )
        AUEs[it] = numpy.abs(c - e).mean()
        taus[it], p = scipy.stats.kendalltau( e, c)
        Rs[it], p = scipy.stats.pearsonr( e, c)
        maxEs[it] = max(abs( e-c))
        pers[it] = CorrectSign(c, e)
        r2s[it] = correl(e, c)

    return [avgerr, avgerrs.std()], [RMS, RMSerrs.std()], [AUE, AUEs.std()],[tau, taus.std()], [R, Rs.std()], [maxE, maxEs.std()], [per, pers.std()], [r2, r2s.std()]
     

# Metho for making plots comparing calculated and experimental Data
def ComparePlot(x, y, Title, XLabel, YLabel, xerr, yerr, labels, fileName = 'compare.pdf', limits = None, leg = [1.02, 0.98, 2, 1],expError = 1.0, wOption = 1, symbols = ['ro','bs','gD','rx','bo','gs']): 
    """ Input:
        x, y, xerr, yerr = list of arrays to be plotted 
        Title, XLabel, YLabel = strings for labeling plot
        labels = list of labels for the legend
        fileName = string, optional, save plot to file
        limits = list with [lowlimit, highlimit] default used if None
        leg = list of legend settings, in order, both bbox_to_anchor placement, location number, and number of columns 
        expError = shaded region of plot
    """
    rcParams.update(JCAMDdict(wOption))

    # If limits not provided find the low and high
    if limits == None:
        # adjustment so limits are aleast as big as error bars 
        maxXerr = np.max([i for k in xerr for i in k])
        maxYerr = np.max([i for k in yerr for i in k])
        adjust = np.max([maxXerr, maxYerr]) + 0.5
        # Find limits find the minimum and maximum of all data
        allvals = [i for k in x for i in k]
        allvals += [i for k in y for i in k]
        lowLim = min(allvals) - adjust
        highLim = max(allvals) + adjust
    # Otherwise use the ones provided
    else:
        lowLim = limits[0]
        highLim = limits[1]

    # Set up so that multiple sets of data can be passed, so if only one set is passed make it a list of lists still
    try:
        len(x[0])
    except:
        x = [x]
        y = [y]
        xerr = [xerr]
        yerr = [yerr]
        
    # Initial Plot settings
    fig1 = plt.figure(1)
    ylim(lowLim, highLim)
    xlim(lowLim, highLim)
    plt.xlabel(XLabel)
    plt.ylabel(YLabel)
    plt.title(Title)
    ax1 = fig1.add_subplot(111)
    handles = []

    # Add data to plot for each set
    for i in range(len(x)):
        p1 = ax1.errorbar(x[i],y[i], xerr = xerr[i], yerr = yerr[i], fmt = symbols[i], label = labels[i], capsize = 0.5)
        handles.append(p1)

    # Insert range for "correct" prediction
    L = numpy.array([lowLim-20, highLim +20])
    ax1.plot(L,L, 'k-', markersize = 6.0)
    ax1.fill_between(L, L+expError, L-expError, facecolor = 'wheat', alpha = 0.5, label = "Agree within %.1f" % expError)
    yellow = patches.Patch(color = 'wheat', label = r'$\pm$ %0.1f log units' % expError)
    handles.append(yellow)

    # Add legend
    ax1.legend(bbox_to_anchor = (1.02, 0.98), loc = 2, ncol = 1, borderaxespad = 0., handles = handles)

    savefig(fileName)

    plt.close(fig1)

# ===============================================================================
# Methods from uncertain_check.py David L. Mobley wrote for the SAMPL4 analysis
# ===============================================================================

def normal( y):
    """Return unit normal distribution value at specified location."""
    return 1./sqrt(2*pi) * exp( -y**2/2. )

def compute_range_table( stepsize = 0.001, maxextent = 10  ):
    """Compute integrals of the unit normal distribution and return these tabulated. Returns:
- range: NumPy array giving integration range (x) where integration range runs -x to +x
- integral: NumPy arrange giving integrals over specified integration range.
Arguments (optional):
- stepsize: Step size to advance integration range by each trial. Default: 0.001
- maxextent: Maximum extent of integration range
"""
    #Calculate integration range
    x = arange( 0, maxextent, stepsize )  #Symmetric, so no need to do negative values.
    
    #Calculate distribution at specified x values
    distrib = normal(x)

    integral = zeros( len(x), float)
    for idx in range(1, len(x)):
        integral[idx] = 2*scipy.integrate.trapz( distrib[0:idx+1], x[0:idx+1] ) #Factor of 2 handles symmetry

    return x, integral 

def get_range( integral, rangetable, integraltable):
    """Use rangetable and integral table provided (i.e. from compute_range_table) to find the smallest range of integration for which the integral is greater than the specified value (integral). Return this range as a float."""

    idx = where( integraltable > integral)[0]
    return rangetable[ idx[0]]


#[DLM]Precompute integral of normal distribution so I can look up integration range which gives desired integral
#integral_range, integral = compute_range_table()



def fracfound_vs_error( calc, expt, dcalc, dexpt, integral_range, integral):
    """
    Takes in calculated and experimental values, their uncertainties as well as 
    """
    #Fraction of Gaussian distribution we want to compute
    X = arange( 0, 1.0, 0.01) 
    Y = zeros( len(X))

    for (i, x) in enumerate(X):
        #Determine integration range which gives us this much probability
        rng = get_range( x, integral_range, integral)
        #print x, rng       
 
        #Loop over samples and compute fraction of measurements found
        y = 0.
        #for n in range(0, len(DGcalc)):
        #    sigma_eff = sqrt( sigma_calc[n]**2 + sigma_expt[n]**2 )
        #    absdiff = abs( DGcalc[n] - DGexpt[n])
        #    #print absdiff, n, sigma_eff
        #    if absdiff < rng * sigma_eff: #If the difference falls within the specified range of sigma values, then this is within the range we're looking at; track it
        #        #print "Incrementing y for n=%s, x = %.2f" % (n, x)
        #        y += 1./len(DGcalc)
        #Rewrite for speed
        sigma_eff = sqrt( array(dcalc)**2 + array(dexpt)**2)
        absdiff = sqrt( (array(calc) - array(expt))**2)
        idx = where( absdiff < rng*sigma_eff)[0] 
        Y[i] = len(idx) * 1./len(calc)

    #print Y
    #raw_input()

    return X, Y

