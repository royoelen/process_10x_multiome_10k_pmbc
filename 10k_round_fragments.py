"""
10k_round_fragments.py

This script is use to round atac fragments from CellRanger ARC, to represent fragment instead of read counts
authors: Martijn van der Werf, Roy Oelen

"""

# import libraries
import gzip
import sys

# set up path
infile = f'./10k_PBMC_Multiome_nextgem_Chromium_X_atac_fragments.tsv.gz'
outfile = f'./10k_PBMC_Multiome_nextgem_Chromium_X_rounded_atac_fragments.tsv.gz'
# setup filehandle
fh = gzip.open(outfile,'wt')
with gzip.open(infile, 'rt') as infile:
    while True:  
        # read the line
        line = infile.readline()
        # write metadata
        if line.startswith('#'):
            fh.write(line)
        # split line
        else:
            line = line.strip().split('\t')
            # round so they are counted fragments instead of reads, meaning pairs are counted instead as a single entry
            if int(line[4]) %2:
                line[4] = (int(line[4]) + 1) / 2
            # write to the new file
            fh.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+str(int(line[4]))+"\t"+"\n")
        if not line:
            break

# close filehandle
fh.close()
