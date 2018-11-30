from yambopy import *

def modifyTimeReversalString(fname):
    with open(fname) as f:
        lines = []
        for l in f:
            if l.startswith('#RmTimeRev'):
                print 'removed # from RmTimeRev field'
                lines.append(l[1:])
            else:
                lines.append(l)
    return lines

def writeLines(fname,lines):
    f = open(fname,'w')
    for l in lines:
        f.write(l)
        f.close

def removeTimeReversal(fname):
    lines = modifyTimeReversalString(fname)
    writeLines(fname,lines)

def fixSymm(kfold,fieldDirection):
    """
    Create the FixSymm folder (in the kfold path but only if it is not already
    present) with the new SAVE folder compatible with the symmetries broken by
    the fieldDirection perturbation
    """
    fixSymmFold = kfold+'/FixSymm'
    if not os.path.isdir(fixSymmFold):
        # create the ypp input file
        yppname = kfold+'ypp.in'
        if os.path.isfile(yppname):
            osStr = "rm %s"%yppname
            print "remove file : %s"%yppname
        osStr = "cd %s; ypp -y"%kfold
        os.system(osStr)
        y = YamboIn(filename=yppname)
        y['Efield1'] = fieldDirection
        y.write(yppname)
        # run ypp again to correctly format the file
        os.system(osStr)
        # remove time reversal
        removeTimeReversal(yppname)
        # execute ypp to create the new SAVE folder
        osStr = "cd %s; ypp -F ypp.in"%kfold
        os.system(osStr)
        # execute yambo without options (setup procedure)
        osStr = "cd %s; yambo"%fixSymmFold
        os.system(osStr)
    else:
        print 'FixSymm folder already created'
    return fixSymmFold

# new routines here.....



def makeTDinput(fold,fname,fieldDirection,fieldInt,fieldFreq,fieldWidth,RTstep,NETime,RTbands,RT_CPU):
    y = YamboIn('yambo_rt -q p -v ip -V qp',folder=fold)
    """
    Build the input file for a TD simulation in the "independt particle" approximation.
    Set the relevant parameters for the field and the simulation options.
    """
    # field paramters
    y['Field1_kind'] = 'QSSIN'
    y['Field1_pol'] = 'linear'
    y['Field1_Dir'] = fieldDirection
    y['Field1_Int'] = [fieldInt,'kWLm2']
    y['Field1_Freq'] = [[fieldFreq,fieldFreq],'eV']
    y['Field1_Width'] = [fieldWidth,'fs']
    # simulation parameters
    y['IOtime'] = [[1.0,5.0,0.5],'fs'] #(J,P,CARRIERs - GF - OUTPUT)
    y['RTstep'] = [RTstep,'as']
    y['NETime'] = [NETime,'fs']
    y['RTBands'] = RTbands
    # parameters from yambopy tutorial rt_si.py
    y['GfnQP_Wv']   = [0.05,0.00,0.00]    # Constant damping valence
    y['GfnQP_Wc']   = [0.05,0.00,0.00]    # Constant damping conduction
    y['GfnQP_E']    = [0.00, 1.00, 1.00]  # [EXTQP BSK BSS] E parameters  (c/v) eV|adim|adim

    # parallelization parameters
    y['RT_CPU'] = RT_CPU
    y.write(fold+'/'+fname)

def runYambo_rt(folder,filename,jobname,mpi,omp):
    """
    Run a single Yambo_rt computation and delete the jobname folder is exsists
    """
    jobDirPath = folder+'/'+jobname
    if os.path.isdir(jobDirPath):
        print 'delete '+ jobDirPath
        os.system("rm -r %s"%jobDirPath)
    osString = "cd %s ; "%folder
    osString += "OMP_NUM_THREADS=%d mpirun -np %d yambo_rt -F %s -J %s -C %s"%(omp,mpi,filename,jobname,jobname)
    print 'execute : '+osString
    os.system(osString)
    print 'done!'

def makeYPP_rtOccupationInput(fold,RTbands,yppTimeStep):
    """"
    fold : the folder in which the input file is created
    This function execute ypp_rt -n o e the ypp.in. Then load the input file in y with
    y = YamboIn(filename=fname), modifies the parameters
    and execute ypp -s b -V qp again to produce the correct final file (called ypp.in)
    """
    # if ypp.in exists is removed
    if os.path.isfile(fold+'/ypp.in'):
        osStr = "rm %s/ypp.in"%fold
        print "remove file : %s/ypp.in"%fold
        os.system(osStr)
    # run ypp -s b -V qp to build the input file
    osStr = "cd %s; ypp_rt -n o e -V qp"%fold
    print osStr
    os.system(osStr)
    fname = fold+'/ypp.in'
    y = YamboIn(filename=fname)
    y['QPkrange'][0][3:5] = RTbands
    y['TimeStep'][0] = 10.0
    #print y
    y.write(fname)
    # run again ypp_rt
    os.system(osStr)

def runYPP_rt(fold,filename,jobname,outfold):
    """
    Run a single YPP_rt computation and delete the outfold folder is exsists
    jobname : the name of the folder with the results of yambo_rt used as input
    """
    jobDirPath = fold+'/'+outfold
    if os.path.isdir(jobDirPath):
        print 'delete '+ jobDirPath
        os.system("rm -r %s"%jobDirPath)
    osString = "cd %s; OMP_NUM_THREADS=1 mpirun -np 1 ypp_rt -F %s -J %s -C %s"%(fold,filename,jobname,outfold)
    print 'execute : '+osString
    os.system(osString)
    print 'done!'
