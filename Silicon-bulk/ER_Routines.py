import shutil as sh
import os

def mkDirFromFilePath(filepath,outDir):
    """"
    Create the folder that contain the file specified with its filepath.
    If the folder already exists do not do anything
    """
    pathList = filepath.split(os.sep)
    folder = outDir
    for level in pathList[0:len(pathList)-1]:
        folder += os.sep+level
        if not os.path.isdir(folder):
            print 'create folder : ', folder
            os.mkdir(folder)
    return folder

def copySAVEfolder(SAVEpath,outDir):
    """"
    Copy the SAVE folder in the outDir path (only if it not already exists)
    """
    newfoldPath = outDir+os.sep+SAVEpath
    if not os.path.isdir(newfoldPath):
        print 'copy ' ,SAVEpath , 'in ', newfoldPath
        sh.copytree(SAVEpath,newfoldPath)

def writeTerm(s,write):
    if write:
        print s

def createResultsDirectory(searchFolder,resultDir,skipFoundFiles = True, write = False):
    if not os.path.isdir(resultDir):
        os.mkdir(resultDir)

    for sf in searchFolder:
        print ''
        print 'search in folder : ', sf
        for fold, subdirs, files in os.walk(sf):
            for f in files:
                filepath = fold + os.sep

                # extract the 'o-*', '.log' output files and 'r_setup' files
                if f.find('o-') != -1 or f.endswith('log') or f.find('r_setup') != -1:
                    outpath = mkDirFromFilePath(filepath,resultDir)
                    if skipFoundFiles:
                        if not os.path.isfile(outpath+os.sep+f):
                            #print 'copy ',fold+os.sep+f, ' in ', outpath
                            writeTerm('copy '+fold+os.sep+f+' in '+outpath,write)
                            sh.copy(filepath+f,outpath)
                        else:
                            #print 'file ',outpath+os.sep+f, ' already found!'
                            writeTerm('file '+outpath+os.sep+f+' already found!',write)
                    else:
                        #print 'copy ',fold+os.sep+f, ' in ', outpath
                        writeTerm('copy '+fold+os.sep+f+' in '+outpath,write)
                        sh.copy(filepath+f,outpath)

            #extract the SAVE folder in the yambo path
            for d in subdirs:
                dirpath = fold + os.sep + d
                if dirpath.find('SAVE') != -1 and dirpath.find('yambo') != -1:
                    copySAVEfolder(dirpath,resultDir)
