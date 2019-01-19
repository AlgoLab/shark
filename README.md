# Read-filter
Tool for mapping-free classification of reads, using Bloom filter.

## Dependencies
Read-filter requires the following libraries and tools:

- [sdsl-lite v2.1.1](https://github.com/simongog/sdsl-lite/tree/v2.1.1)
- [KMC v3.1.0](https://github.com/refresh-bio/KMC/tree/v3.1.0)
- [htslib v1.9](https://github.com/samtools/htslib/tree/1.9)

This repository comes with them as submodules.

## Download and Installation
To install the tool, run the following steps.

First, clone the repository and move into it.
```shell
git clone --recursive https://github.com/tmunimib/read-filter.git
cd read-filter
```

If you have KMC3, sdsl-lite, and htslib already installed you can skip the following commands.
```shell
chmod +x Library_installer.sh
./Library_installer.sh
```

Now you can compile read-filter from the root of you local copy of the repository.
```shell
cd <path-to-read-filter-local-repo>
make
```

## Usage
```
read-filter [-v] -t <transcripts> -1 <sample1> [-2 <sample2>] [-k <kmer_size>] [-c <confidence>]
      -h, --help            display this help and exit
      -t, --transcripts     transcripts in FASTA format (can be gzipped)
      -1, --sample1         sample in FASTA/Q (can be gzipped)
      -2, --sample2         second sample in FASTA/Q (can be gzipped)
      -k, --kmer-size       size of the kmers to index (default:31)
      -c, --confidence      confidence for associating a read to a gene (default: 20)
      -b, --bf-size         bloom filter size in GB (default:1)
      -v, --verbose         verbose mode
```