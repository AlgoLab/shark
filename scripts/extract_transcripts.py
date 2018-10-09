import sys, os, gzip

from Bio import SeqIO

def main():
    fa_path = sys.argv[1]
    chrom = sys.argv[2]

    with gzip.open(fa_path, "rt") as fa_handle:
        for record in SeqIO.parse(fa_handle, "fasta"):
            if 'chromosome:GRCh38:{}:'.format(chrom) in record.description:
                SeqIO.write(record, sys.stdout, "fasta")

if __name__ == '__main__':
    main()
