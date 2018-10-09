import sys, os

from Bio import SeqIO

def main():
    fa_path = sys.argv[1]
    min_len = int(sys.argv[2])

    total_len = 0

    for record in SeqIO.parse(fa_path, "fasta"):
        seq = record.seq
        l = len(seq)
        if l < min_len:
            continue
        total_len += l
        descr = record.description.split(" ")
        transcript = descr[0]
        gene = descr[3]
        if gene.startswith('gene'):
            gene = gene.split(':')[1]
        else:
            gene = "unk_gene"
            for elem in descr:
                if elem.startswith('gene'):
                    gene = elem.split(':')[1] 
        record.description = ""
        record.id = "{}_{}".format(transcript, gene)
        SeqIO.write(record, sys.stdout, "fasta")
    print("Total len: {}".format(total_len), file=sys.stderr)
        

if __name__ == '__main__':
    main()
