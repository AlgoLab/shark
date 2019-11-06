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
      -1, --sample1                     sample in FASTA/Q (can be gzipped)
      -2, --sample2                     second sample in FASTA/Q (optional, can be gzipped)
      -o, --out1                        first output sample in FASTA/Q (default: sharked_sample.1)
      -p, --out2                        second output sample in FASTA/Q (default: sharked_sample.2)
      -k, --kmer-size                   size of the kmers to index (default:17, max:31)
      -c, --confidence                  confidence for associating a read to a gene (default:0.6)
      -b, --bf-size                     bloom filter size in GB (default:1)
      -q, --min-base-quality            minimum base quality (assume FASTQ Illumina 1.8+ Phred scale, default:0, i.e., no filtering)
      -s, --single                      report an association only if a single gene is found
      -t, --threads                     number of threads (default:1)
      -v, --verbose                     verbose mode
```

## Example
```
cd examples
tar xvfz chrY.tar.gz
../shark -r chrY/chrY.cdna.fa -1 chrY/chrY_1.fq -2 chrY/chrY_2.fq -o chrY_1.sharked.fq -p chrY_2.sharked.fq > chrY.ssv
```
