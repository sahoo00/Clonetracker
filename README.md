# Clonetracker
Analysis of barcodes of clones of cells in clonetracker experiment

## General Information

fbp1 sequence CCGACCACCGAACGCAACGCACGCA is used to locate the clone barcode ID in the reverse strand. Unique Cell barcode ID is extracted from first 28 sequences of the forward strand. The Cell barcode ID and Clone barcode ID is combined and stored in a fastq file. The fastq files is filtered using quality information. First 16 sequences from this fastq file is the Cell barcode ID and barcode ID is last 48 sequences. The Cell barcodes and Clone barcodes are clustered using DNACLUST tool. A summary of the clusters of the Cell barcodes and Clone barcodes is collected using following table:

Columns in the final output:

CellID - representative sequence ID of a cell

NumReads = Number of reads supporting a cell cluster

NumClones = Number of clones found inside one cell cluster

BigCloneSize = Size of biggest clone inside a cell

Ratio = ratio of NumReads by NumClones

CellBC = 10x Barcode of the cell

CloneID = ID of the clone


1. Extracting barcodes:
```
python clonetrack.py
```

2. Filtering barcodes:
```
python filter.py
```

3. Extracting cells:
```
python convert.py
```

4. Clustering cells:

```
#!/bin/bash
clust=/booleanfs/sahoo/softwares/SequencingData/dnaclust_linux_release3/dnaclust

#$clust GP1/cells-2.fasta -s 0.99 -k 5 > GP1/cells-2.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' GP1/cells-2.cls > GP1/cells-2.stats
#$clust GP1/barcode-2.fasta -s 0.99 -k 5 > GP1/barcode-2.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' GP1/barcode-2.cls > GP1/barcode-2.stats
#
#$clust GP1/cells-3.fasta -s 0.99 -k 5 > GP1/cells-3.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' GP1/cells-3.cls > GP1/cells-3.stats
#$clust GP1/barcode-3.fasta -s 0.99 -k 5 > GP1/barcode-3.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' GP1/barcode-3.cls > GP1/barcode-3.stats
#
#$clust GP1/cells-4.fasta -s 0.99 -k 5 > GP1/cells-4.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' GP1/cells-4.cls > GP1/cells-4.stats
#$clust GP1/barcode-4.fasta -s 0.99 -k 5 > GP1/barcode-4.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' GP1/barcode-4.cls > GP1/barcode-4.stats

#$clust PE1/cells-2.fasta -s 0.99 -k 5 > PE1/cells-2.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' PE1/cells-2.cls > PE1/cells-2.stats
#$clust PE1/barcode-2.fasta -s 0.99 -k 5 > PE1/barcode-2.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' PE1/barcode-2.cls > PE1/barcode-2.stats
#
#$clust PE1/cells-3.fasta -s 0.99 -k 5 > PE1/cells-3.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' PE1/cells-3.cls > PE1/cells-3.stats
#$clust PE1/barcode-3.fasta -s 0.99 -k 5 > PE1/barcode-3.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' PE1/barcode-3.cls > PE1/barcode-3.stats

#$clust PE1/cells-4.fasta -s 0.99 -k 5 > PE1/cells-4.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' PE1/cells-4.cls > PE1/cells-4.stats
#$clust PE1/barcode-4.fasta -s 0.99 -k 5 > PE1/barcode-4.cls
#awk 'BEGIN{OFS=FS="\t"}{print $1,NF}' PE1/barcode-4.cls > PE1/barcode-4.stats
```

5. Collecting results:
```
python collect.py
```

References:
https://manuals.cellecta.com/clonetracker-xp-lentiviral-barcode-libraries/

