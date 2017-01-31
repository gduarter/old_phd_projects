#!/usr/bin/env python
import pickle, sys
import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import matplotlib.patches as mpatches # for error ranges
from scipy import stats # for linear regressions
from sklearn.metrics import mean_squared_error
from matplotlib import gridspec
import tools as t

with open('data.txt','w+') as f:
    freeSolv = pickle.load(open('database.pickle'))

    dG = []
    dGerr = []
    expt = []
    dexpt = []
    dH = []
    dHerr = []
    dS = []
    dSerr = []

    for key in freeSolv.keys():
        dG.append(freeSolv[key]['calc'])            #kJ/mol
        dGerr.append(freeSolv[key]['d_calc'])
        expt.append(freeSolv[key]['expt'])
        dexpt.append(freeSolv[key]['expt'])
        dH.append(freeSolv[key]['calc_h'])
        dHerr.append(freeSolv[key]['d_calc_h'])
        dS.append(freeSolv[key]['calc_s (cal/mol.K)'])          #kJ/K.mol
        dSerr.append(freeSolv[key]['d_calc_s (cal/mol.K)'])

    dG = np.array(dG) / 0.239006
    dGerr = np.array(dGerr) / 0.239006
    expt = np.array(expt) / 0.239006
    dexpt = np.array(dexpt) / 0.239006
    dH = np.array(dH) / 0.239006
    dHerr = np.array(dHerr) / 0.239006
    dS = np.array(dS) / 0.239006
    dSerr = np.array(dSerr) / 0.239006

    lim1 = -100
    lim2= 20
#    # Plot DG_experimental X DG_FreeSolv
#    pl.rc('font', family='serif')
#    pl.rc('xtick', labelsize='x-small')
#    pl.rc('ytick', labelsize='x-small')
#    slope_new, intercept_new, r_new, p_new, std_new = stats.linregress(expt,dG)
#    rms = np.sqrt(mean_squared_error(expt,dG))
#    fig = pl.figure(figsize=(3.25,3.25))
#    ax = fig.add_subplot(111)
#    ax.set_ylim(lim1,lim2)
#    ax.set_xlim(lim1,lim2)
#    ax.plot([-20/0.239006,20/0.239006],[-20/0.239006,20/0.239006], lw=3)
#    points = ax.errorbar(expt, dG, yerr=dGerr, fmt="o", markersize=4)
#    #legend1 = pl.legend(points,[r"RMSE = %.3f $kcal\cdot mol^{-1}$; $R^2$ = %.3f" % (rms,r_new**2)],loc="upper left")
#    #ax.legend(loc="upper left",prop={"size":20})
#    ax.set_xlabel(r"$\Delta G^{hyd}_{expt}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
#    ax.set_ylabel(r"$\Delta G^{hyd}_{FreeSolv}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
#    patch = mpatches.Patch(color='blue', alpha=0.5, label=r'$\mathrm{\mathsf{within}}$ $\pm\,5.0\,kJ \cdot mol^{-1}$')
#    ax.fill_between([lim1,lim2],[lim1+5,lim2+5],[lim1-5,lim2-5],facecolor = 'blue', alpha = 0.5, label = r'$\pm$ 1.0 $kcal\cdot mol^{-1}')
#    ax.legend( handles=[patch],loc='upper left',prop={'size':10})
#    pl.xticks(fontsize=8)
#    pl.yticks(fontsize=8)
#    #pl.gca().add_artist(legend1)
#    pl.tight_layout()
#    fig.savefig("experimental_X_freeSolv.pdf")
#    ###
#    data = t.stats_array( dG, expt, np.zeros(len(expt)), 100, noise=False )
#    f.write("dG_expt x dG_calc plot\n")
#    f.write("Average error = %.4f +- %.4f\n" % (data[0][0], data[0][1]))
#    f.write("RMS = %.4f +- %.4f\n" % (data[1][0], data[1][1]))
#    f.write("AUE = %.4f +- %.4f\n" % (data[2][0], data[2][1]))
#    f.write("Kendall tau = %.4f +- %.4f\n" % (data[3][0], data[3][1]))
#    f.write("Pearson R = %.4f +- %.4f\n" % (data[4][0], data[4][1]))
#    f.write("\n")

    # Plot DG_FreeSolv X DH
    pl.rc('font', family='serif')
    pl.rc('xtick', labelsize=14)
    pl.rc('ytick', labelsize='xx-small')
    slope_new, intercept_new, r_new, p_new, std_new = stats.linregress(dH,dG)
    fig = pl.figure(figsize=(7,4))
    gs = gridspec.GridSpec(1,5)
    ax1 = fig.add_subplot(gs[0,0:2])#121)
    ax1.set_title("(a)",loc='left',fontdict={'fontsize': 8,
                 'fontweight' : pl.rcParams['axes.titleweight'],
                 'verticalalignment': 'baseline',
                 'horizontalalignment': 'left'})
    ax1.set_ylim(lim1-7,lim2+7)
    #pl.gca().set_aspect('equal')
    points = ax1.errorbar(dH, dG, xerr=dHerr, yerr=dGerr, fmt="o", markersize=4)
    #legend1 = pl.legend(points,[r"$R^2$ = %.3f" % r_new**2],loc="upper left",)
    #ax.legend(loc="upper left",prop={"size":8})
    ax1.set_xlabel(r"$\Delta H^{hyd}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
    ax1.set_ylabel(r"$\Delta G_{FreeSolv}^{hyd}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
    pl.xticks(np.arange(-160, 0, 30),fontsize=8)
    pl.yticks(fontsize=8)
    #pl.gca().add_artist(legend1)
    # Plot DG_FreeSolv X DS
    slope_new, intercept_new, r_new, p_new, std_new = stats.linregress(dS,dG)
    ax2 = fig.add_subplot(gs[0,2:],sharey=ax1)#122)
    ax2.set_title("(b)", loc='left',fontdict={'fontsize': 8,
    'fontweight' : pl.rcParams['axes.titleweight'],
    'verticalalignment': 'baseline',
    'horizontalalignment': 'left'})
    #ax.set_ylim(lim1,lim2)
    #pl.gca().set_aspect('equal')
    #pl.axis('scaled')
    points = ax2.errorbar(dS, dG, xerr=dSerr, yerr=dGerr, fmt="o", markersize=4)
    #legend2 = pl.legend(points,[r"$R^2$ = %.3f" % r_new**2],loc="upper left",)
    #ax.legend(loc="upper left",prop={"size":8})
    ax2.set_xlabel(r"$\Delta S^{hyd}$ ($J\cdot mol^{-1} \cdot K^{-1}$)", fontsize=10)
    #ax2.set_ylabel(r"$\Delta G_{FreeSolv}^{hyd}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
    pl.xticks(fontsize=8)
    pl.yticks(fontsize=8)
    #pl.gca().add_artist(legend2)
    # Minimize overlap and save
    pl.tight_layout()
    fig.savefig("enthalpy_and_entropy.pdf")

#    ###
#    data = t.stats_array( dH, dG, np.zeros(len(expt)), 100, noise=False )
#    f.write("dH x dG plot\n")
#    f.write("Average error = %.4f +- %.4f\n" % (data[0][0], data[0][1]))
#    f.write("RMS = %.4f +- %.4f\n" % (data[1][0], data[1][1]))
#    f.write("AUE = %.4f +- %.4f\n" % (data[2][0], data[2][1]))
#    f.write("Kendall tau = %.4f +- %.4f\n" % (data[3][0], data[3][1]))
#    f.write("Pearson R = %.4f +- %.4f\n" % (data[4][0], data[4][1]))
#    f.write("\n")
#    ###
#    data = t.stats_array( dS, dG, np.zeros(len(expt)), 100, noise=False )
#    f.write("dS x dG plot\n")
#    f.write("Average error = %.4f +- %.4f\n" % (data[0][0], data[0][1]))
#    f.write("RMS = %.4f +- %.4f\n" % (data[1][0], data[1][1]))
#    f.write("AUE = %.4f +- %.4f\n" % (data[2][0], data[2][1]))
#    f.write("Kendall tau = %.4f +- %.4f\n" % (data[3][0], data[3][1]))
#    f.write("Pearson R = %.4f +- %.4f\n" % (data[4][0], data[4][1]))
#    f.write("\n")

#    # Plot bar plot with conformational enthalpy, solvation enthalpy of 10 molecules
#    # Biggest conformational enthalpies
#    # Pre-plot definitions
#    data = [[freeSolv[key]['h_conf'],key] for key in freeSolv.keys()]
#    plotSample = sorted(data, key=lambda l:l[0],reverse=True)[:10]
#    sampleData = [freeSolv[plotSample[i][1]] for i in range(len(plotSample))]
#    nSamples = len(sampleData)
#    wantedEntries = {'h_conf': r'$\Delta H_{hyd}^{conf}$',
#                 'calc_s (cal/mol.K)': r'$T/Delta S_{hyd}'}
#    nSamples = len(sampleData)*len(wantedEntries.keys())
#    widthBar = 250
#    widthSpace = 0.1
#    colorList = ['b','r']
#    indSpace = 1
#    step = indSpace/2.
#    pos = np.array([30, 63, 96, 129, 162, 195, 228, 261, 294,
#                    327, 360, 393, 426, 459, 492, 525, 558,
#                    591, 624, 657, 690])/0.7
#    # Plot
#    fig = pl.figure()
#    pl.axis("off")
#    pl.title("(a)",loc="left")
#    ax = fig.add_axes([0.15, 0.15, 0.65, 0.7])
#    ax.spines['right'].set_color('none')
#    ax.spines['top'].set_color('none')
#    ax.xaxis.set_ticks_position('bottom')
#    ax.yaxis.set_ticks_position('right')
#    ax.tick_params(axis='y', direction='out', width=0, length=0,
#              labelsize=16, pad=0, color='white')
#    ax.tick_params(axis='x', labelsize=14)
#    ax.spines['left'].set_linewidth(3)
#    ax.spines['bottom'].set_linewidth(3)
#    ax.spines['left'].set_position('zero')
#    ax.spines['bottom'].set_position('zero')
#    j = 0
#    y_ticks = []
#    while j < nSamples:
#        for elem in sampleData:
#            for i, entry in enumerate(wantedEntries):
#                if entry == 'calc_s (cal/mol.K)':
#                    ax.barh(pos[j]/0.239006+i*20, elem[entry] *298.15/1000./0.239006, widthBar,
#                            edgecolor='k', linewidth=1,
#                            facecolor=colorList[i])
#                    y_ticks.append(pos[j]+i*20)
#                elif entry == 'h_conf':
#                    ax.barh(pos[j]/0.239006+i*20, elem[entry]/0.239006, widthBar,
#                            edgecolor='k', linewidth=1,
#                            facecolor=colorList[i])
#                    y_ticks.append(pos[j]+i*20)
#                else:
#                    ax.barh(pos[j]/0.239006+i*20, elem[entry]/0.239006, widthBar,
#                            edgecolor='k', linewidth=1,
#                            facecolor=colorList[i])
#                    y_ticks.append(pos[j]+i*20)
#                if i == 1: pos = pos + 100 # remember list indices
#                j += 1
#    # Find tick positions
#    pairs = np.array(y_ticks).reshape((len(y_ticks)/2,2))
#    y_set = []
#    for elem in pairs:
#        y_set.append((elem[0]+elem[1])/2./0.239006)
#    # Start plot
#    ax.set_yticks(y_set)
#    ax.set_yticklabels([r'$\leftarrow$'+freeSolv[plotSample[i][1]]['nickname']
#                        for i in range(len(plotSample))], size='medium',minor=False)
#    red = mpatches.Patch(fc='red',ec='black',label=r'$\Delta H^{hyd}_{conf}$')
#    #green = mpatches.Patch(fc='green',ec='black',label=r'$\Delta G_{hyd}^{vdW}$')
#    blue = mpatches.Patch(fc='blue',ec='black',label=r'$T\Delta S^{hyd}$')
#    #orange = mpatches.Patch(fc='orange',ec='black',label=r'$\Delta G_{hyd}^{q}$')
#    handles = [red,blue]
#    matplotlib.rcParams['legend.frameon'] = False
#    #matplotlib.rcParams['figure.figsize'] = [12.0,6.0]
#    ax.legend(bbox_to_anchor=(0.5,1.15), loc='upper center', ncol=2, borderaxespad=0.4,
#          handles=handles, fontsize=20)
#    ax.set_xlabel(r'$kJ\cdot mol^{-1}$',size=20,labelpad=10)
#    ax.xaxis.set_ticks(np.arange(-80, 70, 20.0))
#    pl.gca().set_aspect(0.026)
#    pl.savefig('barh_positiveHconf.svg',transparent=True)
#    pl.savefig('barh_positiveHconf.pdf')
#
#    ###########################
#    # Plot Enthalpy calc X expt
#    common = pickle.load(open('common.pickle'))
#    use = np.genfromtxt('experimental_enthalpies.txt',usecols=0,dtype=str)
#    common = pickle.load(open('common_update.pickle'))
#    dH_fs = []
#    dHerr_fs = []
#    dH_or = []
#    dHerr_or = []
#    dG_fs = []
#    dGerr_fs = []
#    dG_ex = []
#    for elem in freeSolv:
#        for key in common.keys():
#            if key == freeSolv[elem]['iupac']:
#                common[key]['calc_H'] = freeSolv[elem]['calc_h']
#                common[key]['d_calc_H'] = freeSolv[elem]['d_calc_h']
#                common[key]['FreeSolvEntry'] = elem
#                common[key]['calc_G'] = freeSolv[elem]['calc']
#                common[key]['d_calc_G'] = freeSolv[elem]['d_calc']
#                common[key]['expt'] = freeSolv[elem]['expt']
#            else: pass
#    #pickle.dump(common, open('common_update.pickle','wb'))
#    for key in common.keys():
#        if key in use:
#            dH_fs.append(common[key]['calc_H'])
#            dHerr_fs.append(common[key]['d_calc_H'])
#            dH_or.append(common[key]['expt_H'])
#            dHerr_or.append(common[key]['d_expt_H'])
#            dG_fs.append(common[key]['calc_G'])
#            dGerr_fs.append(common[key]['d_calc_G'])
#            dG_ex.append(common[key]['expt'])
#        else:
#            continue
#    # Convert to array
#    dH_fs = np.array(dH_fs)/0.239006
#    dHerr_fs = np.array(dHerr_fs)/0.239006
#    dH_or = np.array(dH_or)
#    dHerr_or = np.array(dHerr_or)
#    dG_fs = np.array(dG_fs)/0.239006
#    dGerr_fs = np.array(dGerr_fs)/0.239006
#    dG_ex = np.array(dG_ex)/0.239006
#    # Plot DH_expt x DH_FreeSolv
#    pl.rc('font', family='serif')
#    meanError = np.mean(np.absolute(dH_fs - dH_or))
#    slope_new, intercept_new, r_new, p_new, std_new = stats.linregress(dH_fs,dH_or)
#    rms = np.sqrt(mean_squared_error(dH_fs,dH_or))
#    fig = pl.figure(figsize=(6.7,3.0))
#    ax = fig.add_subplot(121)
#    points = ax.errorbar(dH_or, dH_fs, xerr=dHerr_or, yerr=dHerr_fs, fmt="bo", markersize=4)
#    ax.plot([-80,0],[-80,0],"k-",lw=2)
#    ax.set_xlim([-80,0])
#    ax.set_ylim([-80,0])
#    start, end = ax.get_xlim()
#    ax.xaxis.set_ticks(np.arange(start, end, 10.0))
#    patch = mpatches.Patch(color='blue', alpha=0.5, label=r'within $\pm\,1.0\,kcal\cdot mol^{-1}$')
#    ax.fill_between([start,end],[start+1/0.239006,end+1/0.239006],[start-1/0.239006,end-1/0.239006],facecolor = 'blue', alpha = 0.5, label = r'$\pm$ 1.0 $kcal\cdot mol^{-1}')
#    annotation = r"RMSE = %.3f $kJ\cdot mol^{-1}$;" % rms
#    annotation += "\n"
#    annotation +=r"$R^2$ = %.3f;" % r_new**2
#    annotation += "\n"
#    annotation += "mean error = %.3f $kJ\cdot mol^{-1}$" % meanError
#    #legend1 = pl.legend(points,[annotation],loc="lower right")
#    #ax.legend(loc="upper left",prop={"size":8})
#    ax.set_xlabel(r"$\Delta H^{hyd}_{expt}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
#    ax.set_ylabel(r"$\Delta H^{hyd}_{FreeSolv}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
#    ax.set_title("(a)",loc='left',fontdict={'fontsize': 8,
#                 'fontweight' : pl.rcParams['axes.titleweight'],
#                 'verticalalignment': 'baseline',
#                 'horizontalalignment': 'left'})
#    pl.xticks(fontsize=8)
#    pl.yticks(fontsize=8)
#    #pl.gca().add_artist(legend1)
#    # Plot DG_expt x DG_FreeSolv
#    meanError2 = np.mean(np.absolute(dG_fs - dG_ex))
#    slope_new2, intercept_new2, r_new2, p_new2, std_new2 = stats.linregress(dG_fs,dG_ex)
#    rms2 = np.sqrt(mean_squared_error(dG_fs,dG_ex))
#    ax2 = fig.add_subplot(122)
#    points2 = ax2.errorbar(dG_ex, dG_fs, yerr=dGerr_fs, fmt="bo", markersize=4)
#    ax2.plot([-25,20],[-25,20],"k-",lw=2)
#    ax2.set_xlim([-25,20])
#    ax2.set_ylim([-25,20])
#    start, end = ax2.get_xlim()
#    ax2.xaxis.set_ticks(np.arange(start, end, 5.0))
#    patch = mpatches.Patch(color='blue', alpha=0.5, label=r'$\mathrm{mathsf{within}}$ $\pm\,1.0\,kcal\cdot mol^{-1}$')
#    ax2.fill_between([start,end],[start+1/0.239006,end+1/0.239006],[start-1/0.239006,end-1/0.239006],facecolor = 'blue', alpha = 0.5, label = r'$\pm$ 1.0 $kcal\cdot mol^{-1}')
#    annotation2 = r"RMSE = %.3f $kJ\cdot mol^{-1}$;" % rms2
#    annotation2 += "\n"
#    annotation2 +=r"$R^2$ = %.3f;" % r_new2**2
#    annotation2 += "\n"
#    annotation2 += "mean error = %.3f $kJ\cdot mol^{-1}$" % meanError2
#    #legend2 = pl.legend(points2,[annotation2],loc="lower right")
#    #ax2.legend(loc="upper left",prop={"size":8})
#    ax2.set_xlabel(r"$\Delta G^{hyd}_{expt}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
#    ax2.set_ylabel(r"$\Delta G^{hyd}_{FreeSolv}$ ($kJ\cdot mol^{-1}$)", fontsize=10)
#    ax2.set_title("(b)",loc='left',fontdict={'fontsize': 8,
#                  'fontweight' : pl.rcParams['axes.titleweight'],
#                  'verticalalignment': 'baseline',
#                  'horizontalalignment': 'left'})
#    pl.xticks(fontsize=8)
#    pl.yticks(fontsize=8)
#    #ax2.legend( handles=[patch],loc='lower right',prop={'size':8})
#    #pl.gca().add_artist(legend2)
#    # Minimize overlap and save
#    pl.tight_layout()
#    fig.savefig("FreeSolv_and_ORCHYD.pdf")
#    ####
#    data = t.stats_array( dH_fs, dH_or, np.zeros(len(expt)), 100, noise=False )
#    f.write("dH_ex x dH_fs plot\n")
#    f.write("Average error = %.4f +- %.4f\n" % (data[0][0], data[0][1]))
#    f.write("RMS = %.4f +- %.4f\n" % (data[1][0], data[1][1]))
#    f.write("AUE = %.4f +- %.4f\n" % (data[2][0], data[2][1]))
#    f.write("Kendall tau = %.4f +- %.4f\n" % (data[3][0], data[3][1]))
#    f.write("Pearson R = %.4f +- %.4f\n" % (data[4][0], data[4][1]))
#    f.write("\n")
#    data = t.stats_array( dG_fs, dG_ex, np.zeros(len(expt)), 100, noise=False )
#    f.write("dG_exp x dG_calc plot\n")
#    f.write("Average error = %.4f +- %.4f\n" % (data[0][0], data[0][1]))
#    f.write("RMS = %.4f +- %.4f\n" % (data[1][0], data[1][1]))
#    f.write("AUE = %.4f +- %.4f\n" % (data[2][0], data[2][1]))
#    f.write("Kendall tau = %.4f +- %.4f\n" % (data[3][0], data[3][1]))
#    f.write("Pearson R = %.4f +- %.4f\n" % (data[4][0], data[4][1]))
#    f.write("\n")
#









