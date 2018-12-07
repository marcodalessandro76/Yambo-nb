from yambopy import *
import GW_Routines as GW
from copy import deepcopy

def readLinesFromFile(fname):
    with open(fname) as f:
        lines = []
        for l in f:
            lines.append(l)
    return lines

def writeLines(fname,lines):
    f = open(fname,'w')
    for l in lines:
        f.write(l)
        f.close

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

def removeTimeReversal(fname):
    lines = modifyTimeReversalString(fname)
    writeLines(fname,lines)

def addStringToFile(fname,string):
    lines = readLinesFromFile(fname)
    lines.append(string)
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

def make_rt_input_file(path,fname='rt_default.in',**kwargs):
    """
    Build the input file for a RT simulation in the independent particle approximation.
    """
    y = YamboIn('yambo_rt -q p -v ip -V qp',folder=path)
    for k,v in kwargs.iteritems():
        y[k] = v
    y.write(path+'/'+fname)

def build_rt_dictionary(dic,path,rt_par):
    """
    Update the yambo dictionary with the parameters of the choosen RT computations and make the associated
    RT input file. The outputFile name does not contain the .carriers or .external_field part since all the
    output files can be managed in this way.
    """
    fieldInt = dic.keys()
    for f in fieldInt:
        dic_rt = deepcopy(rt_par)
        dic_rt['Field1_Int'] = [f,'kWLm2']
        radical = 'rt_int'+str(int(f))+'_freq'+str(rt_par['Field1_Freq'][0][0])+'_step'+str(rt_par['RTstep'][0])
        inputFile = radical+'.in'
        jobName = radical
        outputFile = path+'/'+jobName+'/'+'o-'+radical
        dic[f]['rt'] = {'inpf' : inputFile, 'jobn' : jobName, 'outf' : outputFile, 'parameters' : dic_rt}

        make_rt_input_file(path,fname = inputFile, **dic_rt)

def run_yambo_rt(path,filename,jobname,mpi,omp):
    """
    Run a single Yambo_rt computation and delete the jobname folder is exsists
    """
    jobDirPath = path+'/'+jobname
    if os.path.isdir(jobDirPath):
        print 'delete '+ jobDirPath
        os.system("rm -r %s"%jobDirPath)
    osString = "cd %s ; "%path
    osString += "OMP_NUM_THREADS=%d mpirun -np %d yambo_rt -F %s -J %s -C %s"%(omp,mpi,filename,jobname,jobname)
    print 'execute : ',osString

    os.system(osString)
    print 'done!'

def run_rt(dic,path,mpi,omp,skip = False):
    """
    Run a bunch of RT simulations
    """
    for f,y in dic.iteritems():
        carriersOut = y['rt']['outf']+'.carriers'
        if skip:
            if os.path.isfile(carriersOut):
                print 'skip the computation for : ',carriersOut
            else:
                run_yambo_rt(path,y['rt']['inpf'],y['rt']['jobn'],mpi,omp)
        else:
            run_yambo_rt(folder,y['rt']['inpf'],y['jobn'],mpi,omp)

def parser_rt_carriers(fname):
    """
    Read the carriers outputfile of the RT simulation
    """
    print 'parsing file : ',fname
    larray = GW.parserArrayFromFile(fname)
    time = []
    dnElec = []
    dnHoles = []
    for l in larray:
        time.append(l[0])
        dnElec.append(l[2])
        dnHoles.append(l[3])

    return {'time':time,'dnElec':dnElec,'dnHoles':dnHoles}

def parser_rt_external_field(fname):
    """
    Read the eternal_field outputfile of the RT simulation
    """
    print 'parsing file : ',fname
    larray = GW.parserArrayFromFile(fname)
    time = []
    eInt = []
    fluence = []
    for l in larray:
        time.append(l[0])
        eInt.append(l[7])
        fluence.append(l[8])

    return {'time':time,'field_int':eInt,'fluence':fluence}

def get_rt_results(dic):
    """
    Update the yambo dictionary with the results of the RT simulation
    """
    fieldInt = dic.keys()
    for f,y in dic.iteritems():
        carriers_outf = y['rt']['outf']+'.carriers'
        carriers_res = parser_rt_carriers(carriers_outf)
        ext_field_outf = y['rt']['outf']+'.external_field'
        ext_field_res = parser_rt_external_field(ext_field_outf)
        y['rt']['results'] = {'time':carriers_res['time'],
                              'field_int':ext_field_res['field_int'],
                              'fluence':ext_field_res['fluence'],
                              'dnElec':carriers_res['dnElec'],
                              'dnHoles':carriers_res['dnHoles']
                             }

def make_ypp_rt_noe_input(path,rt_bands,ypp_time_step):
    """"
    This function execute ypp_rt -n o and generate the file ypp.in in the path folder. Then load the input
    file in y with y = YamboIn(filename=fname), modifies the parameters and execute ypp -s b -V qp again to
    produce the correct final file
    """
    # if ypp.in exists is removed
    if os.path.isfile(path+'/ypp.in'):
        osStr = "rm %s/ypp.in"%path
        print "remove file : %s/ypp.in"%path
        os.system(osStr)
    # run ypp -s b -V qp to build the input file
    osStr = "cd %s; ypp_rt -n o e -V qp"%path
    print osStr
    os.system(osStr)
    fname = path+'/ypp.in'
    y = YamboIn(filename=path+'/ypp.in')
    y['QPkrange'][0][3:5] = rt_bands
    y['TimeStep'][0] = ypp_time_step
    y.write(fname)
    # run again ypp_rt
    os.system(osStr)

def build_ypp_noe(path,dic,ypp_time_step=10.0):
    make_ypp_rt_noe_input(path,dic['parameters']['RTBands'],ypp_time_step)

    radical = 'int'+str(int(dic['parameters']['Field1_Int'][0]))+'_freq'+str(dic['parameters']['Field1_Freq'][0][0])+'_step'+str(dic['parameters']['RTstep'][0])
    inputFile = 'ypp_noe_'+radical+'.in'
    outputFile = path+'/'+'ypp_noe_'+radical+'/'+'o-rt_'+radical+'.YPP-RT_occupations_DATA'
    dic['noe'] = {'jobn':'ypp_noe_'+radical,'inpf':inputFile,'outf':outputFile}

def run_ypp_rt(path,filename,jobname,outfold):
    """
    Run a single YPP_rt computation in the path folder and delete the outfold folder is exsists
    jobname : the name of the folder with the results of yambo_rt used as input
    """
    outDirPath = path+'/'+outfold
    if os.path.isdir(outDirPath):
        print 'delete '+ outDirPath
        os.system("rm -r %s"%outDirPath)
    osString = "cd %s; OMP_NUM_THREADS=1 mpirun -np 1 ypp_rt -F %s -J %s -C %s"%(path,filename,jobname,outfold)
    print 'execute : '+osString
    os.system(osString)
    print 'done!'

def run_noe(path,dic):
    jobname = dic['jobn']
    outfold = dic['noe']['jobn']
    run_ypp_rt(path,'ypp.in',jobname,outfold)

def parser_ypp_noe(fname):
    """
    Read the ypp_noe outputfile and return the energy and the occupations
    at the latest time
    """
    print 'parsing file : ',fname
    larray = GW.parserArrayFromFile(fname)
    energy = []
    occupations = []
    for l in larray:
        energy.append(l[0])
        occupations.append(l[-1])

    return {'energy':energy,'occupations':occupations}

def get_ypp_noe_results(dic):
    results = parser_ypp_noe(dic['noe']['outf'])
    dic['noe']['results'] = results

def make_cs_rt_input(path,fname = 'cs_default.in',**kwargs):
    y = YamboIn('yambo_rt -b -k hartee -g n -p c -V qp',folder=path)
    for k,v in kwargs.iteritems():
        y[k] = v
    y.write(path+'/'+fname)

def build_cs_rt(dic,path,par,time = 100):
    """
    Update the yambo dictionary with the key associated to the cs computation.Add to the par
    dictionary the reference to the neq carriers and make the input file
    """
    for f,y in dic.iteritems():
        rt_par = y['rt']['parameters']
        radical  = 'cs_int'+str(int(f))+'_freq'+str(rt_par['Field1_Freq'][0][0])+'_time'+str(time)
        inputFile = radical+'.in'
        jobName = radical
        outputFile = path+'/'+radical+'/'+'o-'+radical+'.qp'
        y['cs'] = {}
        y['cs'][time] = {'inpf' : inputFile, 'jobn' : jobName, 'outf' : outputFile}

        par['GfnRTdb'] = "f @ %s fs < %s/ndb.RT_carriers"%(time,y['rt']['jobn'])
        par['XfnRTdb'] = "f @ %s fs < %s/ndb.RT_carriers"%(time,y['rt']['jobn'])

        make_cs_rt_input(path,fname = inputFile,**par)

def run_cs_rt(dic,path,mpi,omp,skip = False):
    """
    Run a bunch of cs RT simulations
    """
    for f,y in dic.iteritems():
        val = y['cs'].values()[0] #extract the dictionary associated to the time key
        if skip:
            if os.path.isfile(val['outf']):
                print 'skip the computation for : ',val['outf']
            else:
                run_yambo_rt(path,val['inpf'],val['jobn'],mpi,omp)
        else:
            run_yambo_rt(folder,y['cs']['inpf'],val['jobn'],mpi,omp)

def get_cs_rt_results(dic):
    """
    Reads the output of the cs RT calculations and add the appropriate fields in the yambo
    dictionary
    """
    for f,y in dic.iteritems():
        val = y['cs'].values()[0] #extract the dictionary associated to the time key
        print 'read file : ', val['outf']
        y['cs']['results'] = GW.parserCOHSEXout(val['outf'])
