from qepy import *
from yambopy import *
from qepy import PwIn

def SiInputFile():
    """
    Define a Quantum espresso base input file for silicon
    """
    qe = PwIn()
    qe.atoms = [['Si',[0.125,0.125,0.125]],
                ['Si',[-.125,-.125,-.125]]]
    qe.atypes = {'Si': [28.086,"Si.pbe-mt_fhi.UPF"]} # pbe pseudo
    #qe.atypes = {'Si': [28.086,"Si.vbc.UPF"]} # vbc pseudo

    qe.control['wf_collect'] = '.true.'
    qe.control['pseudo_dir'] = "'./pseudos'"
    qe.system['celldm(1)'] = 10.3
    qe.system['occupations'] = "'fixed'"
    qe.system['nat'] = 2
    qe.system['ntyp'] = 1
    qe.system['ibrav'] = 2
    qe.electrons['conv_thr'] = 1e-8
    return qe

def scfSimulation(scfDict,k,e):
    """
    Manage the dictionary and write the input file for a single scf
    simulation. To be called only from the buildScf method
    """
    qe = SiInputFile()
    qe.control['calculation'] = "'scf'"
    qe.kpoints = [k,k,k]
    qe.system['ecutwfc'] = e
    fileName = 'k'+str(k)+'_ecut'+str(e)
    qe.control['prefix'] = "'scf/output/"+fileName+"'"
    qe.write('scf/input/'+fileName+'.scf')

    scfDict[k][e] = {
        'inputFile' : 'scf/input/'+fileName+'.scf',
        'outputFile' : 'scf/output/'+fileName+'.log'}

def buildScf(kpoints,ecut):
    """
    Build the QE directory structure and input files for a bunch of scf computations.
    Build the scf dictionary
    """
    scfDict = {}
    if not os.path.isdir('scf'):
        os.mkdir('scf')
    if not os.path.isdir('scf/input'):
        os.mkdir('scf/input')
    if not os.path.isdir('scf/output'):
        os.mkdir('scf/output')
    for k in kpoints:
        scfDict[k] = {}
        for e in ecut:
            scfSimulation(scfDict,k,e)

    return scfDict

def runPw(inputFile,outputFile,nthreads):
    """
    Run a single Pw simulation.
    """
    pw = 'pw.x'
    runString =  "mpirun -np %d %s -inp %s > %s"%(nthreads,pw,inputFile,outputFile)
    print 'execute : '+runString
    os.system(runString)

def runScf(dic,nthreads,skip = False):
    """
    Run a bunch of scf simulations, one for each value of kpoints and ecut
    """
    kval = dic.keys()
    kval.sort()
    for k in kval:
        evalues = dic[k].keys()
        evalues.sort()
        for e in evalues:
            if skip:
                if os.path.isfile(dic[k][e]['outputFile']):
                    print 'skip the computation for : '+dic[k][e]['outputFile']
                else:
                    runPw(dic[k][e]['inputFile'],dic[k][e]['outputFile'],nthreads)
            else:
                runPw(dic[k][e]['inputFile'],dic[k][e]['outputFile'],nthreads)


# Parser of the QE output file to extract the total energy
def get_energyStringFromFile(file):
    with open(file, 'r') as f:
        for line in f:
            if '!    total energy' in line:
                energyString = line
    return energyString

def extractEnergy(energyString):
    lenTotal = len(energyString)
    start = len('!    total energy              =     ')
    stop = lenTotal - len(' Ry')
    energy = float(energyString[start:stop])
    return energy

def get_totalEnergy(file):
    energyString = get_energyStringFromFile(file)
    energy = extractEnergy(energyString)
    return energy

def nscfSimulation(dic,k,e,n):
    """
    Manage the dictionary and write the input file for a single nscf
    simulation. To be called only from the buildNscf method
    """
    qe = SiInputFile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.system['force_symmorphic'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = n
    qe.kpoints = [k,k,k]
    qe.system['ecutwfc'] = e
    qe.control['prefix'] = "'nscf/output/k%s_ecut%s_nb%s'"%(k,e,n)

    fileName = 'k'+str(k)+'_ecut'+str(e)+'_nbnd'+str(n)
    qe.write('nscf/input/'+fileName+'.nscf')

    dic[k][e][n] = {
        'inputFile' : 'nscf/input/'+fileName+'.nscf',
        'outputFile' : 'nscf/output/'+fileName+'.log',
        'outFolder' : 'nscf/output/k%s_ecut%s_nb%s'%(k,e,n)+'.save'}

def buildNscf(kpoints,ecut,nb,kconv,ecutconv):
    """
    Build the QE directory structure and input files for a bunch of nscf computations.
    Build the nscf dictionary.
    """
    dic = {}
    if not os.path.isdir('nscf'):
        os.mkdir('nscf')
    if not os.path.isdir('nscf/input'):
        os.mkdir('nscf/input')
    if not os.path.isdir('nscf/output'):
        os.mkdir('nscf/output')
    for k in kpoints:
        dic[k] = {}
        for e in ecut:
            dic[k][e] = {}
            for n in nb:
                dic[k][e][n] = {}
                #If the nscf.save folder is missing copy the scf.save folder associated to the converged scf computation and
                #rename so to have the same name of the nscf run
                if not os.path.isdir("nscf/output/k%s_ecut%s_nb%s.save"%(k,e,n)):
                    cpString = "cp -r scf/output/k%s_ecut%s.save nscf/output/"%(kconv,ecutconv)
                    print 'execute : '+cpString
                    os.system(cpString)
                    cpString = "mv  nscf/output/k%s_ecut%s.save nscf/output/k%s_ecut%s_nb%s.save"%(kconv,ecutconv,k,e,n)
                    print 'execute : '+cpString
                    os.system(cpString)

                nscfSimulation(dic,k,e,n)

    return dic

def runNscf(dic,nthreads,skip = False):
        """
        Run a bunch of nscf simulation, one for each value of kpoints, ecut and nbnds
        """
        kval = dic.keys()
        kval.sort()
        for k in kval:
            evalues = dic[k].keys()
            evalues.sort()
            for e in evalues:
                nb = dic[k][e].keys()
                nb.sort()
                for n in nb:
                    if skip:
                        if os.path.isfile(dic[k][e][n]['outputFile']):
                            print 'skip the computation for : '+dic[k][e][n]['outputFile']
                        else:
                            runPw(dic[k][e][n]['inputFile'],dic[k][e][n]['outputFile'],nthreads)
                    else:
                        runPw(dic[k][e][n]['inputFile'],dic[k][e][n]['outputFile'],nthreads)

def runP2y(dic):
    """
    Perform p2y and execute yambo without options in al the .save folders of the nscf simulations
    """
    kval = dic.keys()
    kval.sort()
    for k in kval:
        evalues = dic[k].keys()
        evalues.sort()
        for e in evalues:
            nb = dic[k][e].keys()
            nb.sort()
            for n in nb:
                osString = "cd %s;p2y;yambo"%dic[k][e][n]['outFolder']
                print 'execute : '+osString
                os.system(osString)


def nscfOutFolderSplit(val):
    out = val.partition('nscf/output/')[2]
    out = out.partition('.save')[0]
    return out

def buildYambo(dic):
    """
    Take the nscfDict as input and build the yambo directory structure and
    the base yambo dictionary for a bunch of yambo computations.
    """
    if not os.path.isdir('yambo'):
        os.mkdir('yambo')
    yamboDic = {}
    for k in dic:
        yamboDic[k] = {}
        for e in dic[k]:
            yamboDic[k][e] = {}
            for n in dic[k][e]:
                folderName = nscfOutFolderSplit(dic[k][e][n]['outFolder'])
                if not os.path.isdir('yambo/'+folderName):
                    os.mkdir('yambo/'+folderName)
                    #copy the SAVE folder
                    osString = "cp -r %s/SAVE yambo/%s" %(dic[k][e][n]['outFolder'],folderName)
                    print 'execute : ' + osString
                    os.system(osString)

                #create the HF dic with folder key
                yamboDic[k][e][n] = {'folder' : 'yambo/'+folderName}

    return yamboDic

def buildHFinput(fold,fname,exRL):
    y = YamboIn('yambo -x -V All',folder=fold)
    #y['EXXRLvcs'] = exRL
    #y['EXXRLvcs'] = str(exRL)+' mHa'
    y['EXXRLvcs'] = [exRL,'mHa']
    y.write(fold+'/'+fname)

def buildHF(ydic,exRL):
    """
    Build the input file for a yambo HF computation and update the yambo dictionary
    with the paramters of the choosen HF computations. Note that the inputFile field
    does not include the path of the file which is specified in the folder field. This
    choice is due to the way in which yambo is called.
    """
    kpoints = ydic.keys()
    for k in kpoints:
        ecuts = ydic[k].keys()
        for e in ecuts:
            nb = ydic[k][e].keys()
            for n in nb:
                ydic[k][e][n]['hf'] = {}
                for ex in exRL:
                    jobname = 'hf_exRL'+str(ex)
                    inpfile = 'hf_exRL'+str(ex)+'.in'
                    outfile = 'o-hf_exRL'+str(ex)+'.hf'
                    buildHFinput(ydic[k][e][n]['folder'],inpfile,ex)
                    ydic[k][e][n]['hf'][ex]= {'inputFile':inpfile,
                    'jobName':jobname,
                    'outputFile':ydic[k][e][n]['folder']+'/'+jobname+'/'+outfile}

def runYambo(folder,filename,jobname,nthreads):
    """
    Run a single Yambo computation and delete the jobname folder is exsists
    """
    jobDirPath = folder+'/'+jobname
    if os.path.isdir(jobDirPath):
        print 'delete '+ jobDirPath
        os.system("rm -r %s"%jobDirPath)
    osString = "cd %s ; "%folder
    osString += "mpirun -np %d yambo -F %s -J %s -C %s"%(nthreads,filename,jobname,jobname)
    print 'execute : '+osString
    os.system(osString)

def runHF(ydic,nthreads,skip = False):
    """
    Run a bunch of HF simulations
    """
    kpoints = ydic.keys()
    kpoints.sort()
    for k in kpoints:
        ecuts = ydic[k].keys()
        ecuts.sort()
        for e in ecuts:
            nb = ydic[k][e].keys()
            nb.sort()
            for n in nb:
                folder = ydic[k][e][n]['folder']
                for y in ydic[k][e][n]['hf'].values():
                    if skip:
                        if os.path.isfile(y['outputFile']):
                            print 'skip the computation for : '+y['outputFile']
                        else:
                            runYambo(folder,y['inputFile'],y['jobName'],nthreads)
                    else:
                        runYambo(folder,y['inputFile'],y['jobName'],nthreads)

def parserArrayFromFile(fname):
    """"
    Build a list that contains the lines of fname avoiding the ones that start with #
    """
    lines = []
    with open(fname) as f:
        for l in f:
            if not l.startswith('#'):
                lines.append(l)
    #split each line in a list (of strings)
    larray = [[] for i in range(len(lines))]
    for ind,l in enumerate(lines):
        larray[ind] = l.split()
    #convert the string to double
    for row in range(len(larray)):
        for col in range(len(larray[row])):
            larray[row][col] = float(larray[row][col])
    return larray

def parserHFout(fname):
    """"
    Return a set of list with the output of the .hf file
    """
    larray = parserArrayFromFile(fname)
    KP = []
    BND = []
    E0 = []
    EHF = []
    DFT = []
    HF = []
    for rows,l in enumerate(larray):
        KP.append(l[0])
        BND.append(l[1])
        E0.append(l[2])
        EHF.append(l[3])
        DFT.append(l[4])
        HF.append(l[5])
    return KP,BND,E0,EHF,DFT,HF

def getHFresults(ydic):
    """
    Reads the output of the HF calculations and add the appropriate fields in the yambo
    dictionary
    """
    kpoints = ydic.keys()
    for k in kpoints:
        ecuts = ydic[k].keys()
        for e in ecuts:
            nb = ydic[k][e].keys()
            for n in nb:
                for exRl,y in ydic[k][e][n]['hf'].iteritems():
                    print 'read file : ' + y['outputFile']
                    y['KP'],y['BND'],y['E0'],y['EHF'],y['DFT'],y['HF'] = parserHFout(y['outputFile'])
