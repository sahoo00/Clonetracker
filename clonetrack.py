from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools
import gzip
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re

def Ham_dist(s1, s2):
    """Calculate Hamming distance between 2 sequences."""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def long_check(pattern, text):
    """Naively and understandably generate (Hamming_dist, start_offset)
    for matches with distance 0 or 1"""
    m = len(pattern)
    for i in range(len(text) - m + 1):
        d = Ham_dist(pattern, text[i:i+m])
        if d < 4:
            yield d, i

def getBarcodes(file_f, file_r, ofile):
    handle = open(ofile, "w")
    count = 0
    index = 0
    fbp1 = "CCGACCACCGAACGCAACGCACGCA"
    f_iter = FastqGeneralIterator(gzip.open(file_f,"rt"))
    r_iter = FastqGeneralIterator(gzip.open(file_r,"rt"))
    for (f_id, f_seq, f_q), (r_id, r_seq, r_q) in zip(f_iter,r_iter):
        index += 1
        m = list(long_check(fbp1, r_seq))
        if (len(m) > 0):
            m = sorted(m, key=lambda x:x[0])
            s1 = m[0][1] + len(fbp1)
            if (s1 + 48) < len(r_seq):
                handle.write("@%s\n%s\n+\n%s\n" 
                    % (f_id, f_seq[:28]+r_seq[s1:(s1+48)], f_q[:28]+r_q[s1:(s1+48)]))
            count += 1
        if count > 10:
            #break
            count += 0
        if (index % 1000) == 0:
            print(index, count)
    handle.close()
    print("{} records written to barcode.fastq".format(count))

def T1():
    file_f = 'Barcode Tracking Experiment/10X Barcode Tissue/GP1_fastq/GP1_CKDL220000660-1a-SI_GA_H7_HGC77DSX3_S4_L003_R1_001.fastq.gz'
    file_r = 'Barcode Tracking Experiment/10X Barcode Tissue/GP1_fastq/GP1_CKDL220000660-1a-SI_GA_H7_HGC77DSX3_S4_L003_R2_001.fastq.gz'
    ofile = "GP1/barcode.fastq"
    getBarcodes(file_f, file_r, ofile)

def T2():
    file_f = 'Barcode Tracking Experiment/Amplified Barcode Seq/GP_1/Renal_Xenograft_Barcoded_GP1_S1_L001_R1_001.fastq.gz'
    file_r = 'Barcode Tracking Experiment/Amplified Barcode Seq/GP_1/Renal_Xenograft_Barcoded_GP1_S1_L001_R2_001.fastq.gz'
    ofile = "GP1/barcode-2.fastq"
    getBarcodes(file_f, file_r, ofile)

def T3():
    file_f = 'Barcode Tracking Experiment/10X Barcode Tissue/PE_fastq/PE2_CKDL220000660-1a-SI_GA_H9_HGC77DSX3_S1_L003_R1_001.fastq.gz'
    file_r = 'Barcode Tracking Experiment/10X Barcode Tissue/PE_fastq/PE2_CKDL220000660-1a-SI_GA_H9_HGC77DSX3_S1_L003_R2_001.fastq.gz'
    ofile = "PE1/barcode.fastq"
    getBarcodes(file_f, file_r, ofile)

def T4():
    file_f = 'Barcode Tracking Experiment/Amplified Barcode Seq/PE_2/Renal_Xenograft_Barcoded_PE2_S2_L001_R1_001.fastq.gz'
    file_r = 'Barcode Tracking Experiment/Amplified Barcode Seq/PE_2/Renal_Xenograft_Barcoded_PE2_S2_L001_R2_001.fastq.gz'
    ofile = "PE1/barcode-2.fastq"
    getBarcodes(file_f, file_r, ofile)

#T1()
#T2()
#T4()
#T3()
