# Shark
Fast tool for mapping-free gene separation of reads, using Bloom filter.

## Dependencies
Shark requires the following libraries and tools:
- [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite/tree/v2.1.1)

This repository comes with it as submodule.

## Download and Installation
To install the tool, run the following steps.

First, clone the repository and move into it.
```shell
git clone --recursive https://github.com/AlgoLab/shark.git
cd shark
cd sdsl-lite
./install.sh ..
cd ..
make
```

## Usage
```
shark [-v] -r <references> -1 <sample1> [-2 <sample2>] [-k <kmer size>] [-c <confidence>] [-b <filter size>] [-q <min base quality>] [-s]
      -h, --help                        display this help and exit
      -r, --reference                   reference sequences in FASTA format (can be gzipped)
      -1, --sample1                     sample in FASTQ (can be gzipped)
      -2, --sample2                     second sample in FASTQ (optional, can be gzipped)
      -o, --out1                        first output sample in FASTQ (default: sharked_sample.1)
      -p, --out2                        second output sample in FASTQ (default: sharked_sample.2)
      -k, --kmer-size                   size of the kmers to index (default:17, max:31)
      -c, --confidence                  confidence for associating a read to a gene (default:0.6)
      -b, --bf-size                     bloom filter size in GB (default:1)
      -q, --min-base-quality            minimum base quality (assume FASTQ Illumina 1.8+ Phred scale, default:0, i.e., no filtering)
      -s, --single                      report an association only if a single gene is found
      -t, --threads                     number of threads (default:1)
      -v, --verbose                     verbose mode
```

## Output format

`shark` outputs to `stdout` a ssv file reporting associations between reads and genes.
The first element of each line is the name of the read whereas the second element is the gene identifier.

Reads in the samples that pass the filter step are stored in the files passed as argument to `-o` and `-p`.

## Example

A small example is provided in the example directory.
* `ENSG00000277117.fa` is the gene sequence of gene ENSG00000277117 in FASTA format
* `sample_1.fq` and `sample_2.fq` are a paired-end RNA-Seq sample simulated from genes ENSG00000277117 and ENSG00000275464

To filter out the reads sequenced from gene ENSG00000277117, run `shark` as follows:

```
./shark -r example/ENSG00000277117.fa -1 example/sample_1.fq \
                                      -2 example/sample_2.fq \
                                      -o example/sharked.sample_1.fq \
                                      -p example/sharked.sample_2.fq > example/ENSG00000277117.ssv
```

The results should be equal to: `example/ENSG00000277117.truth.ssv`, `example/sharked.sample_1.truth.fq`, and  `example/sharked.sample_2.truth.fq`.
