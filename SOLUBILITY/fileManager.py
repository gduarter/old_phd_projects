#!/usr/bin/env python
import subprocess as sub
import shutil, os 
import glob

def make_directory(name):
    """ Creates a directory
    VARIABLE        I/O         Type        Description
    ----------------------------------------------------------------------------
    name            input       string      name of the directory
    """
    try:
        os.mkdir(name)
    except OSError:
        pass

def CreateRunningScripts(name,num):
    """ Fixes LAMBDANR on run scripts and changes the Job name in the running script 

    VARIABLE        I/O         Type        Description
    ----------------------------------------------------------------------------
    name            input       string      name of the sh script (before .X.sh)
    num             input       integer     number of files to be created

    """
    sub.call('cp %s.X.sh %s.0.sh' % (name,name), shell=True)
    ideal_number = num
    for number in range(ideal_number):
        if not os.path.isfile("%s.%s.sh" % (name,number)):
            sub.call("cp %s.0.sh %s.%s.sh" % (name, name, number), shell=True)
    return

def FixLambdasAndNames(name,files="*.sh"):
    """ Fixes LAMBDANR on run scripts and changes the Job name in the running script 

    VARIABLE        I/O         Type        Description
    ----------------------------------------------------------------------------
    name            input       string      name to be given to the job
    files           input       string      file extension

    """
    discrepancy = "export LAMBDANR"
    job_name = '#SBATCH --job-name='
    run_files = glob.glob(files)
    for filename in run_files:
        f = open(filename, "r")
        name_lambdas = filename.split(".")[1]
        lines = f.readlines()
        f.close()
        new_lines = []
        for line in lines:
            if discrepancy in line:
                file_lambdas = line.split()[-1].split("=")[-1]
                if file_lambdas != name_lambdas:
                    new_line = line.replace(file_lambdas, name_lambdas)
                    new_lines.append(new_line)
                else:
                    new_lines.append(line)
            elif job_name in line:
                replacement = '#SBATCH --job-name="%s_%s"\n' % (name, name_lambdas)
                new_lines.append(replacement)
            else:
                new_lines.append(line)
        f2 = open(filename, "w")
        f2.writelines(new_lines)
        f2.close()
    return

def CopyMDP(num):
    """ Replicates mdp files for each lambda value.

    VARIABLE        I/O         Type        Description
    ----------------------------------------------------------------------------
    num            input        integer     number of replicas

    You should run the scripts on GreenPlanet using the following command:
    for ((n=0; n<20; n++)); do sbatch run.${n}.sh; done
    """
    files = ['minimize.X.mdp','equil_nvt.X.mdp','equil_npt.X.mdp','equil_npt2.X.mdp','prod.X.mdp']
    for elem in files:
        f = elem
        file  = open(f, 'r')
        s = file.read()
        N = num # Number of different states to be calculated
        for i in range(N):
            of = f.replace('X', str(i))
            oss = s.replace('XXX', str(i))
            with open(of, 'w') as outputfile:
                outputfile.write(oss)
    return

def OrganizeDGFiles(files="prod.*.xvg", result="results.txt", pickle="results.pickle"):
    """ Organize files in a special directory to simplify analysis

    VARIABLE        I/O         Type        Description
    ----------------------------------------------------------------------------
    files           input       string      production stage xvg files
    results         input       string      text file with alchemical analysis results
    pickle          input       string      pickle file with alchemical analysis results

    """
    xvg_files = glob.glob(files)
    make_directory("RESULTS")
    path = os.getcwd()
    for xvg in xvg_files:
        init = os.path.join("%s","%s" % (path,xvg))
        final = os.path.join("%s","RESULTS/%s" % (path,xvg))
        sutil.move(init, final)
    result_files = [result, pickle]
    for result in result_files:
        init = os.path.join("%s","%s" % (path,result))
        final = os.path.join("%s","RESULTS/%s" % (path,result))
        shutil.move(init, final)
    return

def TerminationCheckDG(logfiles,xvgfiles,Num):
    """ Check if a Free Energy job is finished

    VARIABLE        I/O         Type        Description                                         
    ---------------------------------------------------------------------------------------     
    logfiles        input       list        production stage log filenames stored in a list
    xvgfiles        input       list        production stage xvg filenames stored in a list
    Num             input       integer     number of simulations that were run simultaneously
    incomplete      output      dictionary  dictionary with incomplete simulations entries
    
    """
    prodfiles = ["prod.%s.xvg" % i for i in range(Num)]
    incomplete = {}
    store1 = []
    store2 = []
    if len(xvgfiles) < Num:
        loc = os.getcwd()
        message1 = "Simulation in %s has crashed somewhere" % loc
        for xvgfile in xvgfiles:
            if xvgfile not in prodfiles:
                store1.append(xvgfile)
        incomplete[loc] = {tuple(store1),message1}
        for logfile in logfiles:
            var = sub.check_output("grep 'Finished mdrun' %s | wc -l" % logfile, shell=True)
            loc = os.getcwd()
            if var == "0":
                store2.append(logfile)
        message2 = "These simulations were stopped and need to be resumed"
        incomplete[loc] = {tuple(store2),message2}
        return incomplete
    elif len(xvgfiles) == Num:
        for logfile in logfiles:
            var = sub.check_output("grep 'Finished mdrun' %s | wc -l" % logfile, shell=True)
            loc = os.getcwd()
            if var == "0":
                store2.append(logfile)
                message = "These simulations were stopped and need to be resumed"
                incomplete[loc] = {tuple(store2),message}
        return incomplete
    else:
        return incomplete

def TerminationCheck(logfile):
    """ Check if a regular job is finished

    VARIABLE        I/O         Type        Description                                         
    ---------------------------------------------------------------------------------------
    logfile         input       string      production stage log filenames stored in a list
    message         output      string      

    """
    var = sub.check_output("grep 'Finished mdrun' %s | wc -l" % logfile, shell=True)
    if var == "0":
        return "Simulation incomplete. Restart including -cpi ****.cpt  in the mdrun line"
    else:
        return "Simulation ran smoothly"




    















