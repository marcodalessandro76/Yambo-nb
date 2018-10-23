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
    print 'done!'

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

def nscfSimulation(dic,k,e,nb):
    """
    Manage the dictionary and write the input file for a single nscf
    simulation. To be called only from the buildNscf method
    """
    qe = SiInputFile()
    qe.control['calculation'] = "'nscf'"
    qe.electrons['diago_full_acc'] = ".true."
    qe.system['force_symmorphic'] = ".true."
    qe.electrons['conv_thr'] = 1e-8
    qe.system['nbnd'] = nb
    qe.kpoints = [k,k,k]
    qe.system['ecutwfc'] = e
    qe.control['prefix'] = "'nscf/output/k%s_nb%s'"%(k,nb)

    fileName = 'k'+str(k)+'_nb'+str(nb)
    qe.write('nscf/input/'+fileName+'.nscf')

    dic[k] = {
        'inputFile' : 'nscf/input/'+fileName+'.nscf',
        'outputFile' : 'nscf/output/'+fileName+'.log',
        'outFolder' : 'nscf/output/k%s_nb%s'%(k,nb)+'.save'}

def buildNscf(kpoints,nb,kconv,ecutconv):
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
        #If the nscf.save folder is missing copy the scf.save folder associated to the converged scf computation and
        #rename so to have the same name of the nscf run
        if not os.path.isdir("nscf/output/k%s_nb%s.save"%(k,nb)):
            cpString = "cp -r scf/output/k%s_ecut%s.save nscf/output/"%(kconv,ecutconv)
            print 'execute : '+cpString
            os.system(cpString)
            cpString = "mv  nscf/output/k%s_ecut%s.save nscf/output/k%s_nb%s.save"%(kconv,ecutconv,k,nb)
            print 'execute : '+cpString
            os.system(cpString)

        nscfSimulation(dic,k,ecutconv,nb)

    return dic

def runNscf(dic,nthreads,skip = False):
        """
        Run a bunch of nscf simulation, one for each value of kpoints
        """
        kval = dic.keys()
        kval.sort()
        for k in kval:
            if skip:
                if os.path.isfile(dic[k]['outputFile']):
                    print 'skip the computation for : '+dic[k]['outputFile']
                else:
                    runPw(dic[k]['inputFile'],dic[k]['outputFile'],nthreads)
            else:
                runPw(dic[k]['inputFile'],dic[k]['outputFile'],nthreads)

def runP2y(dic):
    """
    Perform p2y and execute yambo without options in al the .save folders of the nscf simulations
    """
    kval = dic.keys()
    kval.sort()
    for k in kval:
        osString = "cd %s;p2y;yambo"%dic[k]['outFolder']
        print 'execute : '+osString
        os.system(osString)

def nscfOutFolderSplit(val):
    out = val.partition('nscf/output/')[2]
    #out = out.partition('.save')[0]
    out = out.partition('_')[0]
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
        folderName = nscfOutFolderSplit(dic[k]['outFolder'])
        if not os.path.isdir('yambo/'+folderName):
            os.mkdir('yambo/'+folderName)
            #copy the SAVE folder
            osString = "cp -r %s/SAVE yambo/%s" %(dic[k]['outFolder'],folderName)
            print 'execute : ' + osString
            os.system(osString)
        else:
            print 'yambo/'+folderName + ' already present'

        #create the yambo dictionary with folder key
        yamboDic[k] = {'folder' : 'yambo/'+folderName}

    return yamboDic

def makeHFinput(fold,fname,exRL,firstbnd,lastbnd):
    #QPkrange : QP generalized K(points) indices. The format is
    #first k-point|last-kpoint|first-band|last-band|
    #In the dictionary is saved as
    #y[QPkrange] = [[firstk,lastk,firstbnd,lastbnd],'']
    y = YamboIn('yambo -x -V rl',folder=fold)
    y['EXXRLvcs'] = [1000.0*exRL,'mHa']
    krange = y['QPkrange'][0][:2]
    kbandrange = krange + [firstbnd,lastbnd]
    y['QPkrange'] = [kbandrange,'']
    y.write(fold+'/'+fname)

def buildHF(ydic,gcomp,firstbnd,lastbnd):
    """
    Build the input file for a yambo HF computation and update the yambo dictionary
    with the paramters of the choosen HF computations. Note that the inputFile field
    does not include the path of the file which is specified in the folder field. This
    choice is due to the way in which yambo is called.
    """
    kpoints = ydic.keys()
    for k in kpoints:
        ydic[k]['hf'] = {}
        for ex in gcomp:
            jobname = 'hf_gComp'+str(ex)
            inpfile = 'hf_gComp'+str(ex)+'.in'
            outfile = 'o-hf_gComp'+str(ex)+'.hf'
            makeHFinput(ydic[k]['folder'],inpfile,ex,firstbnd,lastbnd)
            ydic[k]['hf'][ex]= {'inputFile':inpfile,
            'jobName':jobname,
            'outputFile':ydic[k]['folder']+'/'+jobname+'/'+outfile}

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
    print 'done!'

def runHF(ydic,nthreads,skip = False):
    """
    Run a bunch of HF simulations
    """
    kpoints = ydic.keys()
    kpoints.sort()
    for k in kpoints:
        folder = ydic[k]['folder']
        for y in ydic[k]['hf'].values():
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
    for rows,l in enumerate(larray):
        KP.append(l[0])
        BND.append(l[1])
        E0.append(l[2])
        EHF.append(l[3])
    return KP,BND,E0,EHF

def getHFresults(ydic):
    """
    Reads the output of the HF calculations and add the appropriate fields in the yambo
    dictionary
    """
    kpoints = ydic.keys()
    for k in kpoints:
        for y in ydic[k]['hf'].values():
            print 'read file : ' + y['outputFile']
            y['KP'],y['BND'],y['E0'],y['EHF'] = parserHFout(y['outputFile'])

def makeCOHSEXinput(fold,fname,gcomp,wg,wn,firstk,lastk,firstbnd,lastbnd):
    #-b static inverse dielectric matrix
    #-k kernel type: hartree
    #-g Dyson Equation solver (n)ewton
    #-p GW approximations (c)OHSEX
    y = YamboIn('yambo -b -k hartee -g n -p c -V qp',folder=fold)
    y['EXXRLvcs'] = [1000.0*gcomp,'mHa']
    y['NGsBlkXs'] = [1000.0*wg,'mHa']
    y['BndsRnXs'] = [1,wn]
    kbandrange = [firstk,lastk] + [firstbnd,lastbnd]
    y['QPkrange'] = [kbandrange,'']
    y.write(fold+'/'+fname)

def buildCOHSEX(ydic,kconv,G0Gconv,wgcomp,wnbnds,firstk,lastk,firstbnd,lastbnd):
    """
    Build the input file for a yambo COHSEX computation and update the yambo dictionary
    with the paramters of the choosen computations. The keys of the ['cs'] dictionary
    are tuple of the form (W_Gcomp,W_nb).
    Note that the inputFile field does not include the path of the file which is specified in
    the folder field. This choice is due to the way in which yambo is called.
    """
    if kconv in ydic.keys():
        ydic[kconv]['cs'] = {}
        for wg in wgcomp:
            for wn in wnbnds:
                jobname = 'cs_wGcomp'+str(wg)+'_wNb'+str(wn)
                inpfile = 'cs_wGcomp'+str(wg)+'_wNb'+str(wn)+'.in'
                outfile = 'o-cs_wGcomp'+str(wg)+'_wNb'+str(wn)+'.qp'
                makeCOHSEXinput(ydic[kconv]['folder'],inpfile,G0Gconv,wg,wn,firstk,lastk,firstbnd,lastbnd)
                ydic[kconv]['cs'][(wg,wn)]= {
                    'inputFile':inpfile,
                    'jobName':jobname,
                    'outputFile':ydic[kconv]['folder']+'/'+jobname+'/'+outfile}
    else:
        print 'k value %s is not present. Add this value to the nscf simulation list'%k

def runCOHSEX(ydic,kconv,nthreads,skip = False):
    """
    Run a bunch of CHOSEX simulations (without empties)
    """
    folder = ydic[kconv]['folder']
    for y in ydic[kconv]['cs'].values():
        if skip:
            if os.path.isfile(y['outputFile']):
                print 'skip the computation for : '+y['outputFile']
            else:
                runYambo(folder,y['inputFile'],y['jobName'],nthreads)
        else:
            runYambo(folder,y['inputFile'],y['jobName'],nthreads)

def parserCOHSEXout(fname):
    """"
    Return a set of list with the output of the .hf file
    """
    larray = parserArrayFromFile(fname)
    KP = []
    BND = []
    E0 = []
    EmE0 = []
    for rows,l in enumerate(larray):
        KP.append(l[0])
        BND.append(l[1])
        E0.append(l[2])
        EmE0.append(l[3])
    return KP,BND,E0,EmE0

def getCOHSEXresults(ydic,kconv):
    """
    Reads the output of the COHSEX calculations and add the appropriate fields in the yambo
    dictionary
    """
    for y in ydic[kconv]['cs'].values():
            print 'read file : ' + y['outputFile']
            y['KP'],y['BND'],y['E0'],y['EmE0'] = parserCOHSEXout(y['outputFile'])

def makePPinput(fold,fname,gcomp,wg,wn,gnbn,firstk,lastk,firstbnd,lastbnd):
    #-d dynamical inverse dielectric matrix
    #-k kernel type: hartree
    #-g Dyson Equation solver (n)ewton
    #-p GW approximations (p)pa
    y = YamboIn('yambo -d -k hartee -g n -p p -V qp',folder=fold)
    y['EXXRLvcs'] = [1000.0*gcomp,'mHa'] #check the name of these parameters
    y['NGsBlkXp'] = [1000.0*wg,'mHa']
    y['BndsRnXp'] = [1,wn]
    y['GbndRnge'] = [1,gnbn]
    kbandrange = [firstk,lastk] + [firstbnd,lastbnd]
    y['QPkrange'] = [kbandrange,'']
    y.write(fold+'/'+fname)

def buildPP(ydic,kconv,G0Gconv,wgconv,wnbnconv,g0nb,firstk,lastk,firstbnd,lastbnd):
    """
    Build the input file for a yambo plasmon pole computation and update the yambo dictionary
    with the paramters of the choosen computations. The keys of the ['pp'] dictionary
    are paramterize the number of empty bands used to compute G0.
    Note that the inputFile field does not include the path of the file which is specified in
    the folder field. This choice is due to the way in which yambo is called.
    """
    if kconv in ydic.keys():
        ydic[kconv]['pp'] = {}
        for gn in g0nb:
            jobname = 'pp_G0nb'+str(gn)
            inpfile = 'pp_G0nb'+str(gn)+'.in'
            outfile = 'o-pp_G0nb'+str(gn)+'.qp'
            makePPinput(ydic[kconv]['folder'],inpfile,G0Gconv,wgconv,wnbnconv,gn,firstk,lastk,firstbnd,lastbnd)
            ydic[kconv]['pp'][gn]= {
                'inputFile':inpfile,
                'jobName':jobname,
                'outputFile':ydic[kconv]['folder']+'/'+jobname+'/'+outfile}
    else:
        print 'k value %s is not present. Add this value to the nscf simulation list'%k

def runPP(ydic,kconv,nthreads,skip = False):
    """
    Run a bunch of plasmon pole simulations
    """
    folder = ydic[kconv]['folder']
    for y in ydic[kconv]['pp'].values():
        if skip:
            if os.path.isfile(y['outputFile']):
                print 'skip the computation for : '+y['outputFile']
            else:
                runYambo(folder,y['inputFile'],y['jobName'],nthreads)
        else:
                runYambo(folder,y['inputFile'],y['jobName'],nthreads)

def getPPresults(ydic,kconv):
    """
    Reads the output of the PP calculations and add the appropriate fields in the yambo
    dictionary
    """
    for y in ydic[kconv]['pp'].values():
            print 'read file : ' + y['outputFile']
            y['KP'],y['BND'],y['E0'],y['EmE0'] = parserCOHSEXout(y['outputFile'])

def makeYPPbandsInput(kfold,study,firstbnd,lastbnd,bands_step,path):
    """"
    kfold : the kpoint folder
    study : the study associated to the chosen ndb.QP file
    This function execute ypp -s b -V qp to produce the ypp.in. Then load the
    input file in y with y = YamboIn(filename=fname), modifies the parameters
    and execute ypp -s b -V qp again to produce the correct final file (called ypp.in)
    """
    # if ypp.in exists is removed
    if os.path.isfile(kfold+'/ypp.in'):
        osStr = "rm %s/ypp.in"%kfold
        print "remove file : %s/ypp.in"%kfold
        os.system(osStr)
    # run ypp -s b -V qp to build the input file
    osStr = "cd %s; ypp -s b -V qp"%kfold
    os.system(osStr)
    fname = kfold+'/ypp.in'
    y = YamboIn(filename=fname)
    y['BANDS_steps'] = bands_step
    kbandrange = y['QPkrange'][0][:2] + [firstbnd,lastbnd]
    y['QPkrange'] = [kbandrange,'']
    if study['jobName'] != 'lda':
        dbname = study['jobName']+'/ndb.QP'
        y['GfnQPdb'] = 'E < '+dbname
    y['BKpts'] = path
    y.write(fname)
    # run again ypp -s b -V qp
    os.system(osStr)

def buildYPPbands(kfold,study,firstbnd,lastbnd,bands_step,path):
    """"
    Build the input file for ypp band calculation to be formed in the kpoints folder
    kfold and associated to the type of computation given by the study variable
    """
    makeYPPbandsInput(kfold,study,firstbnd,lastbnd,bands_step,path)
    sJname = study['jobName']
    study['bnds'] = {'jobName' : sJname+'_bands' ,'outputFile' : 'o-'+sJname+'_bands.bands_interpolated'}

def runYPP(kfold,filename,jobname,nthreads):
    """
    Run a single YPP computation and delete the jobname folder is exsists
    """
    jobDirPath = kfold+'/'+jobname
    if os.path.isdir(jobDirPath):
        print 'delete '+ jobDirPath
        os.system("rm -r %s"%jobDirPath)
    osString = "cd %s; mpirun -np %d ypp -F %s -J %s -C %s"%(kfold,nthreads,filename,jobname,jobname)
    print 'execute : '+osString
    os.system(osString)
    print 'done!'

def runYPPbands(kfold,study,nthreads = 1):
    if 'bnds' in study:
        jobname = study['bnds']['jobName']
        runYPP(kfold,'ypp.in',jobname,nthreads)
    else:
        print "'bnds' key not found"

def parserBandsResults(kfold,study,firstbnd,lastbnd):

    studyfold = study['bnds']['jobName']
    outf = study['bnds']['outputFile']
    fname = kfold+'/'+studyfold+'/'+outf
    print 'parsing file :'+fname
    larray = parserArrayFromFile(fname)

    numbands = lastbnd-firstbnd+1

    kaxis = []
    bndStructure = [[] for i in range(numbands)]
    bndStructure
    for rows,l in enumerate(larray):
        kaxis.append(l[0])
        for ind,b in enumerate(bndStructure):
            b.append(l[ind+1])
    study['bnds']['results'] = {'kaxis' : kaxis, 'bndStructure' : bndStructure}


#######################################################################
def modifyUseBandsString(fname):
    with open(fname) as f:
        lines = []
        for l in f:
            if l.startswith('#UseEbands'):
                print 'removed # from UseBands field'
                lines.append(l[1:])
            else:
                lines.append(l)
    return lines

def writeLines(fname,lines):
    f = open(fname,'w')
    for l in lines:
        f.write(l)
        f.close

def activateUseBand(fname):
    lines = modifyUseBandsString(fname)
    writeLines(fname,lines)
