import pandas as pd
import numpy as np
import re
from collections import Counter
import sys
sys.path.append("/booleanfs2/sahoo/Hegemon/")
sys.path = ["/booleanfs2/sahoo/BoNE/"] + sys.path
import StepMiner as smn
import HegemonUtil as hu
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def getCells(cfile):
    df = pd.read_csv(cfile, header=None, sep="\t")
    df[2] = np.log(df[1])
    df[2].describe()
    thr = hu.getThrData(df[2])
    print(sum(df[2] > thr[0]))
    cells = {k:1 for k in df[df[2] > thr[0]][0]}
    return cells

def getUF(cfile):
    fp = open(cfile, "r")
    uf = {}
    for line in fp:
        line = re.sub("[\r\n]", "", line)
        list1 = line.split("\t")
        for k in list1:
            if k != '':
                uf[k] = list1[0]
    fp.close()
    return uf

def getResults(cfile, ofile, uf, cells):
    idlist = []
    clist = []
    bclist = []
    mlist = []
    fp = open(cfile, "r")
    for line in fp:
        line = re.sub("[\r\n]", "", line)
        list1 = line.split("\t")
        if list1[0] in cells:
            list2 = Counter([uf[k] for k in list1 if k != '' and k in uf])
            idlist += [list1[0]]
            clist += [len(list2)]
            bclist += [len(list1)]
            mlist += [list2.most_common(1)[0][1]]
    fp.close()
    df1 = pd.DataFrame()
    df1['CellID'] = idlist
    df1['NumReads'] = bclist
    df1['NumClones'] = clist
    df1['BigCloneSize'] = mlist
    df1['Ratio'] = df1['NumReads']/df1['NumClones']
    df1.sort_values(by='Ratio')
    df1.to_csv(ofile, sep="\t", index=False)
    return df1

def printCommon(df1, uf):
    chash = {k:1 for k in df1['CellID']}
    list2 = Counter([uf[k] for k in chash if k != '' and k in uf])
    print(list2.most_common(10))

def addCells(file1, ofile, df1, uf):
    handle = open(ofile, "w")
    f_iter = FastqGeneralIterator(open(file1,"rU"))
    index = 0
    df1['CellBC'] = ""
    df1['CloneID'] = ""
    for (f_id, f_seq, f_q) in f_iter:
        id1 = f_id.split(" ")[0]
        if id1 in df1.index:
            index += 1
            df1['CellBC'][id1] = f_seq[:16]
            df1['CloneID'][id1] = uf[id1]
        if (index % 1000) == 0:
            print(index)
    handle.close()
    df1.to_csv(ofile, sep="\t", index=False)

def T1():
    dir1 = "/booleanfs2/sahoo/Data/BooleanLab/Stanford"
    cfile = dir1 + "/GP1/cells-3.stats"
    cells = getCells(cfile)

    dir1 = "/booleanfs2/sahoo/Data/BooleanLab/Stanford"
    cfile = dir1 + "/GP1/barcode-3.cls"
    uf = getUF(cfile)

    dir1 = "/booleanfs2/sahoo/Data/BooleanLab/Stanford"
    cfile = dir1 + "/GP1/cells-3.cls"
    ofile = dir1 + "/GP1/clones-3.txt"
    df1 = getResults(cfile, ofile, uf, cells)
    df1.index = df1['CellID']

    printCommon(df1, uf)
    dir1 = "/booleanfs2/sahoo/Data/BooleanLab/Stanford"
    cfile = dir1 + "/GP1/barcode-3.fastq"
    ofile = dir1 + "/GP1/clones-3.txt"
    addCells(cfile, ofile, df1, uf)

def T2():
    dir1 = "/booleanfs2/sahoo/Data/BooleanLab/Stanford"
    cfile = dir1 + "/PE1/cells-3.stats"
    cells = getCells(cfile)

    dir1 = "/booleanfs2/sahoo/Data/BooleanLab/Stanford"
    cfile = dir1 + "/PE1/barcode-3.cls"
    uf = getUF(cfile)

    dir1 = "/booleanfs2/sahoo/Data/BooleanLab/Stanford"
    cfile = dir1 + "/PE1/cells-3.cls"
    ofile = dir1 + "/PE1/clones-3.txt"
    df1 = getResults(cfile, ofile, uf, cells)
    df1.index = df1['CellID']

    dir1 = "/booleanfs2/sahoo/Data/BooleanLab/Stanford"
    cfile = dir1 + "/PE1/barcode-3.fastq"
    ofile = dir1 + "/PE1/clones-3.txt"
    addCells(cfile, ofile, df1, uf)
    printCommon(df1, uf)

#T1()
T2()
