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
    qe.control['prefix'] = "'nscf/output/k%s_nb%s'"%(k,n)

    fileName = 'k'+str(k)+'_nb'+str(n)
    qe.write('nscf/input/'+fileName+'.nscf')

    dic[k][n] = {
        'inputFile' : 'nscf/input/'+fileName+'.nscf',
        'outputFile' : 'nscf/output/'+fileName+'.log',
        'outFolder' : 'nscf/output/k%s_nb%s'%(k,n)+'.save'}

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
        for n in nb:
            dic[k][n] = {}
            #If the nscf.save folder is missing copy the scf.save folder associated to the converged scf computation and
            #rename so to have the same name of the nscf run
            if not os.path.isdir("nscf/output/k%s_nb%s.save"%(k,n)):
                cpString = "cp -r scf/output/k%s_ecut%s.save nscf/output/"%(kconv,ecutconv)
                print 'execute : '+cpString
                os.system(cpString)
                cpString = "mv  nscf/output/k%s_ecut%s.save nscf/output/k%s_nb%s.save"%(kconv,ecutconv,k,n)
                print 'execute : '+cpString
                os.system(cpString)

            nscfSimulation(dic,k,ecutconv,n)

    return dic


def runNscf(dic,nthreads,skip = False):
        """
        Run a bunch of nscf simulation, one for each value of kpoints and nbnds
        """
        kval = dic.keys()
        kval.sort()
        for k in kval:
            nb = dic[k].keys()
            nb.sort()
            for n in nb:
                if skip:
                    if os.path.isfile(dic[k][n]['outputFile']):
                        print 'skip the computation for : '+dic[k][n]['outputFile']
                    else:
                        runPw(dic[k][n]['inputFile'],dic[k][n]['outputFile'],nthreads)
                else:
                    runPw(dic[k][n]['inputFile'],dic[k][n]['outputFile'],nthreads)

def runP2y(dic):
    """
    Perform p2y and execute yambo without options in al the .save folders of the nscf simulations
    """
    kval = dic.keys()
    kval.sort()
    for k in kval:
        nb = dic[k].keys()
        nb.sort()
        for n in nb:
            osString = "cd %s;p2y;yambo"%dic[k][n]['outFolder']
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
        for n in dic[k]:
            folderName = nscfOutFolderSplit(dic[k][n]['outFolder'])
            if not os.path.isdir('yambo/'+folderName):
                os.mkdir('yambo/'+folderName)
                #copy the SAVE folder
                osString = "cp -r %s/SAVE yambo/%s" %(dic[k][n]['outFolder'],folderName)
                print 'execute : ' + osString
                os.system(osString)
            else:
                print 'yambo/'+folderName + ' already present'

            #create the yambo dictionary with folder key
            yamboDic[k][n] = {'folder' : 'yambo/'+folderName}

    return yamboDic

def buildHFinput(fold,fname,exRL,firstbnd,lastbnd):
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
        nb = ydic[k].keys()
        nb.sort()
        n = nb[0]
        ydic[k][n]['hf'] = {}
        for ex in gcomp:
            jobname = 'hf_gComp'+str(ex)
            inpfile = 'hf_gComp'+str(ex)+'.in'
            outfile = 'o-hf_gComp'+str(ex)+'.hf'
            buildHFinput(ydic[k][n]['folder'],inpfile,ex,firstbnd,lastbnd)
            ydic[k][n]['hf'][ex]= {'inputFile':inpfile,
            'jobName':jobname,
            'outputFile':ydic[k][n]['folder']+'/'+jobname+'/'+outfile}

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
        nb = ydic[k].keys()
        nb.sort()
        n = nb[0]
        folder = ydic[k][n]['folder']
        for y in ydic[k][n]['hf'].values():
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
        nb = ydic[k].keys()
        nb.sort()
        n = nb[0]
        for y in ydic[k][n]['hf'].values():
            print 'read file : ' + y['outputFile']
            y['KP'],y['BND'],y['E0'],y['EHF'] = parserHFout(y['outputFile'])

def buildCOHSEXinput(fold,fname,gcomp,wg,wn,firstk,lastk,firstbnd,lastbnd):
    #-b static inverse dielectric matrix
    #-k kernel type: hartree
    #-g Dyson Equation solver (n)ewton
    #-p GW approximations (c)OHSEX
    y = YamboIn('yambo -b -k hartee -g n -p c -V qp',folder=fold)
    y['EXXRLvcs'] = [1000.0*gcomp,'mHa']
    y['NGsBlkXs'] = [1000.0*wg,'mHa']
    y['BndsRnXs'] = [1,wn]
    #krange = y['QPkrange'][0][:2]
    #kbandrange = krange + [firstbnd,lastbnd]
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
        nb = ydic[kconv].keys()
        nb.sort()
        n = nb[0] #use only the lowest value of nscf_nbnds
        ydic[kconv][n]['cs'] = {}
        for wg in wgcomp:
            for wn in wnbnds:
                jobname = 'cs_wGcomp'+str(wg)+'_wNb'+str(wn)
                inpfile = 'cs_wGcomp'+str(wg)+'_wNb'+str(wn)+'.in'
                outfile = 'o-cs_wGcomp'+str(wg)+'_wNb'+str(wn)+'.qp'
                buildCOHSEXinput(ydic[kconv][n]['folder'],inpfile,G0Gconv,wg,wn,firstk,lastk,firstbnd,lastbnd)
                ydic[kconv][n]['cs'][(wg,wn)]= {
                    'inputFile':inpfile,
                    'jobName':jobname,
                    'outputFile':ydic[kconv][n]['folder']+'/'+jobname+'/'+outfile}
    else:
        print 'k value %s is not present. Add this value to the nscf simulation list'%k

def runCOHSEX(ydic,kconv,nthreads,skip = False):
    """
    Run a bunch of CHOSEX simulations (without empties)
    """
    nb = ydic[kconv].keys()
    nb.sort()
    n = nb[0]
    folder = ydic[kconv][n]['folder']
    for y in ydic[kconv][n]['cs'].values():
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
    nb = ydic[kconv].keys()
    nb.sort()
    n = nb[0]

    for y in ydic[kconv][n]['cs'].values():
            #print ind, y
            print 'read file : ' + y['outputFile']
            y['KP'],y['BND'],y['E0'],y['EmE0'] = parserCOHSEXout(y['outputFile'])

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

def buildCOHSEXWEinput(fold,fname,gcomp,wg,wn,gnbnds,firstbnd,lastbnd):
    #-b static inverse dielectric matrix
    #-k kernel type: hartree
    #-g Dyson Equation solver (n)ewton
    #-p GW approximations (c)OHSEX
    y = YamboIn('yambo -b -k hartee -g n -p c -V all',folder=fold)
    activateUseBand(fold+'/yambo.in')
    y['EXXRLvcs'] = [gcomp,'Ha']
    y['NGsBlkXs'] = [wg,'Ha']
    y['BndsRnXs'] = [1,wn]
    y['GbndRnge'] = [1,gnbnds]
    krange = y['QPkrange'][0][:2]
    kbandrange = krange + [firstbnd,lastbnd]
    y['QPkrange'] = [kbandrange,'']
    #print(y)
    y.write(fold+'/'+fname)

def buildCOHSEXWE(ydic,kconv,G0Gconv,wgconv,wnbndsconv,g0nb,firstbnd,lastbnd):
    """
    Build the input file for a yambo COHSEX (with empties) computation and update the yambo dictionary
    with the paramters of the choosen computations. The keys of the ['cswe'] dictionary
    contain the values of the parameter G0_nb.
    Note that the inputFile field does not include the path of the file which is specified in
    the folder field. This choice is due to the way in which yambo is called.
    """
    if kconv in ydic.keys():
        nb = ydic[kconv].keys()
        nb.sort()
        n = nb[0] #use only the lowest value of nscf_nbnds
        ydic[kconv][n]['cswe'] = {}
        for gn in g0nb:
            jobname = 'cswe_G0nb'+str(gn)
            inpfile = 'cswe_G0nb'+str(gn)+'.in'
            outfile = 'o-cswe_G0nb'+str(gn)+'.qp'
            buildCOHSEXWEinput(ydic[kconv][n]['folder'],inpfile,G0Gconv,wgconv,wnbndsconv,gn,firstbnd,lastbnd)
            ydic[kconv][n]['cswe'][gn]= {
                'inputFile':inpfile,
                'jobName':jobname,
                'outputFile':ydic[kconv][n]['folder']+'/'+jobname+'/'+outfile}
    else:
        print 'k value %s is not present. Add this value to the nscf simulation list'%k

def runCOHSEXWE(ydic,kconv,nthreads,skip = False):
    """
    Run a bunch of CHOSEX simulations (without empties)
    """
    nb = ydic[kconv].keys()
    nb.sort()
    n = nb[0]
    folder = ydic[kconv][n]['folder']
    for y in ydic[kconv][n]['cswe'].values():
        if skip:
            if os.path.isfile(y['outputFile']):
                print 'skip the computation for : '+y['outputFile']
            else:
                runYambo(folder,y['inputFile'],y['jobName'],nthreads)
        else:
            runYambo(folder,y['inputFile'],y['jobName'],nthreads)

def getCOHSEXWEresults(ydic,kconv):
    """
    Reads the output of the COHSEXWE calculations and add the appropriate fields in the yambo
    dictionary
    """
    nb = ydic[kconv].keys()
    nb.sort()
    n = nb[0]

    for y in ydic[kconv][n]['cswe'].values():
            print 'read file : ' + y['outputFile']
            y['KP'],y['BND'],y['E0'],y['EmE0'] = parserCOHSEXout(y['outputFile'])
