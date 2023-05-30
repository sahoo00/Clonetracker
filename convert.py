from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools
import gzip
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def T1(file1, ofile):
    handle = open(ofile, "w")
    f_iter = FastqGeneralIterator(open(file1,"rU"))
    index = 0
    for (f_id, f_seq, f_q) in f_iter:
        index += 1
        handle.write(">%s\n%s\n" % (f_id, f_seq[28:]))
        if (index % 1000) == 0:
            print(index)
    handle.close()

def T2(file1, ofile):
    handle = open(ofile, "w")
    f_iter = FastqGeneralIterator(open(file1,"rU"))
    index = 0
    for (f_id, f_seq, f_q) in f_iter:
        index += 1
        handle.write(">%s\n%s\n" % (f_id, f_seq[:16]))
        if (index % 1000) == 0:
            print(index)
    handle.close()

#T1("GP1/barcode-2.fastq", "GP1/barcode-2.fasta")
#T2("GP1/barcode-2.fastq", "GP1/cells-2.fasta")
#T1("GP1/barcode-3.fastq", "GP1/barcode-3.fasta")
#T2("GP1/barcode-3.fastq", "GP1/cells-3.fasta")
#T1("GP1/barcode-4.fastq", "GP1/barcode-4.fasta")
#T2("GP1/barcode-4.fastq", "GP1/cells-4.fasta")

T1("PE1/barcode-2.fastq", "PE1/barcode-2.fasta")
T2("PE1/barcode-2.fastq", "PE1/cells-2.fasta")
T1("PE1/barcode-3.fastq", "PE1/barcode-3.fasta")
T2("PE1/barcode-3.fastq", "PE1/cells-3.fasta")
#T1("PE1/barcode-4.fastq", "PE1/barcode-4.fasta")
#T2("PE1/barcode-4.fastq", "PE1/cells-4.fasta")

