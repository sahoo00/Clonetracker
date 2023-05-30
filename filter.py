from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
import re

def quality_filter(cfile, ofile):
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    from Bio import SeqIO
    handle = open(ofile, "w")
    for rec in SeqIO.parse(cfile, "fastq"):
        if min(rec.letter_annotations["phred_quality"]) >= 30:
            handle.write(rec.format("fastq"))
    handle.close()


#quality_filter("GP1/barcode-2.fastq", "GP1/barcode-3.fastq")
#quality_filter("GP1/barcode.fastq", "GP1/barcode-4.fastq")
quality_filter("PE1/barcode-2.fastq", "PE1/barcode-3.fastq")
#quality_filter("PE1/barcode.fastq", "PE1/barcode-4.fastq")
