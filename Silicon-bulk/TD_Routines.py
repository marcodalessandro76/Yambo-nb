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
    Create the fixSymm folder (in the kfold path but only if it has not already
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
        print 'fixSymm folder already created'
    return fixSymmFold

def runYambo_rt(folder,filename,jobname,nthreads):
    """
    Run a single Yambo_rt computation and delete the jobname folder is exsists
    """
    jobDirPath = folder+'/'+jobname
    if os.path.isdir(jobDirPath):
        print 'delete '+ jobDirPath
        os.system("rm -r %s"%jobDirPath)
    osString = "cd %s ; "%folder
    osString += "mpirun -np %d yambo_rt -F %s -J %s -C %s"%(nthreads,filename,jobname,jobname)
    print 'execute : '+osString
    os.system(osString)
    print 'done!'
