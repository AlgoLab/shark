# Read-filter
Tool for mapping-free classification of reads, using Bloom filter.

## Dependencies
Read-filter requires the following libraries and tools:
- [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite/tree/v2.1.1)

This repository comes with it as submodule.

## Download and Installation
To install the tool, run the following steps.

First, clone the repository and move into it.
```shell
git clone --recursive https://github.com/AlgoLab/read-filter.git
cd read-filter
cd sdsl-lite/build
./build.sh
cd ../..
make
```

## Usage
```
read-filter [-v] -t <transcripts> -1 <sample1> [-2 <sample2>] [-k <kmer_size>] [-c <confidence>] [-b <bf-size>]
      -h, --help                        display this help and exit
      -r, --reference                   reference sequences in FASTA format (can be gzipped)
      -1, --sample1                     sample in FASTA/Q (can be gzipped)
      -2, --sample2                     second sample in FASTA/Q (optional, can be gzipped)
      -k, --kmer-size                   size of the kmers to index (default:31, max:31)
      -c, --confidence                  confidence for associating a read to a gene (default: 20)
      -b, --bf-size                     bloom filter size in GB (default:1)
      -t, --threads                     number of threads (default:1)
      -v, --verbose                     verbose mode
```

## Example
```
cd examples
tar xvfz chrY.tar.gz
../read-filter -r chrY/chrY.cdna.fa -1 chrY/chrY_1.fq -2 chrY/chrY_2.fq > chrY.ssv
```
