# Scripts to setup exps

##### Prerequisites:
- python3
- biopython
- matplotlib
- wgsim (from samtools package)

### Steps:
1. Get transcripts _Homo_sapiens.GRCh38.cdna.all.fa.gz_ from [here](https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/cdna/)
2. Extract transcripts of a specific chromosome (e.g. 1, 22, Y...)
```bash
python3 extract_transcripts.py [cDNA.fa.gz] [chr] > [chr.fa]
```
3. Extract transcripts longer than [min_len] and fix their headers to simplify downstream analysis, i.e. accuracy computation. This scripts will print the total length of all considered transcripts (to be used in the next step)
```bash
python3 clean_fa.py [chr.fa] [min_len] > [chr_mod.fa]
```
4. Simulate reads from transcripts - [n_reads] is computed by using the formula you can find [here](https://en.wikipedia.org/wiki/Coverage_(genetics)#Calculation) (the length of the genome is the length computed in the previous step):
```bash
wgsim -N [n_reads] -1 [read_len] -e [err_rate] -R 0 -r 0 [chr_mod.fa] [chr.reads.fq] /dev/null
```
5. Run some analysis on coverage:
```bash
grep "^>" [chr_mod.fa] | sed 's/>//g' | while read line ; do grep -c $line [reads.fq] ; done > [chr_transcript_coverage]
python3 draw_histo.py [chr_transcript_coverage]
```
